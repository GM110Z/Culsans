#!/usr/bin/env python3

"""
Resumable targeted KofamScan cell-surface annotation (LPS, O-antigen, and type IV pilus).

Uses calibrated KOfam thresholds, never existing gene/product names.

Input is a text file containing one protein FASTA (.faa) path per genome.
Assembly accessions (GCF_/GCA_) are extracted directly from those paths.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import tarfile

from pathlib import Path
from collections import defaultdict, Counter


ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def save(path, data):
    tmp = path.with_suffix(path.suffix + '.partial')
    tmp.write_text(json.dumps(data, indent=2) + '\n')
    tmp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def prepare(db):
    marker = db / 'prepared_surface.json'

    if marker.exists():
        return json.loads(marker.read_text())

    mapping = defaultdict(list)

    for filename, label in [
        ('lps_biosynthesis_kos.tsv', 'LPS biosynthesis'),
        ('o_antigen_kos.tsv', 'O-antigen sugar biosynthesis'),
        ('pilus_kos.tsv', 'Type IV pilus')
    ]:
        for line in (db / filename).read_text().splitlines():
            pathway, ko = line.split('\t')
            mapping[ko.removeprefix('ko:')].append(label)

    profiles = db / 'surface_profiles'
    profiles.mkdir(exist_ok=True)

    with gzip.open(db / 'ko_list.gz', 'rt') as f:
        text = f.read()

    (db / 'ko_list').write_text(text)

    entries = {
        r['knum']: r
        for r in csv.DictReader(text.splitlines(), delimiter='\t')
    }

    found = set()

    print(
        'Validating archive and extracting pathway profiles',
        flush=True
    )

    # Validate gzip CRC
    with gzip.open(db / 'profiles.tar.gz', 'rb') as f:
        while f.read(4 * 1024 * 1024):
            pass

    with tarfile.open(db / 'profiles.tar.gz', 'r|gz') as archive:
        for member in archive:

            name = Path(member.name).name
            ko = name.removesuffix('.hmm')

            if (
                member.isfile()
                and name.endswith('.hmm')
                and ko in mapping
            ):
                (profiles / name).write_bytes(
                    archive.extractfile(member).read()
                )
                found.add(ko)

    usable = sorted(
        k for k in found
        if entries.get(k, {}).get('threshold', '-') != '-'
    )

    if not usable:
        raise ValueError('No calibrated profiles found')

    hal = db / 'surface.hal'

    hal.write_text(
        ''.join(
            str((profiles / (k + '.hmm')).resolve()) + '\n'
            for k in usable
        )
    )

    info = {
        'mapping': dict(mapping),
        'usable': usable,
        'missing_profiles': sorted(set(mapping) - found),
        'without_threshold': sorted(found - set(usable)),
        'sha256': {
            n: sha(db / n)
            for n in [
                'profiles.tar.gz',
                'ko_list.gz',
                'lps_biosynthesis_kos.tsv',
                'o_antigen_kos.tsv',
                'pilus_kos.tsv'
            ]
        }
    }

    save(marker, info)

    return info


def parse_hits(detail):
    hits = []

    with detail.open() as f:
        for r in csv.reader(f, delimiter='\t'):

            if not r or r[0] != '*':
                continue

            if len(r) != 7:
                raise ValueError(
                    'Unexpected KofamScan output row'
                )

            assembly, protein = r[1].split('|', 1)

            hits.append([
                assembly,
                protein,
                r[2],
                r[3],
                r[4],
                r[5],
                r[6]
            ])

    return hits


def read_faa_list(path):
    """
    Read one FAA path per line and extract assembly accession
    from each path.
    """

    path = path.resolve()

    if not path.is_file():
        raise ValueError(
            f'FAA list does not exist: {path}'
        )

    cohort = []
    inputs = []

    seen_assemblies = set()

    for line in path.read_text().splitlines():

        line = line.strip()

        if not line:
            continue

        faa = Path(line).expanduser()

        # Relative paths are interpreted relative to the
        # directory from which the script is being run.
        if not faa.is_absolute():
            faa = (Path.cwd() / faa).resolve()

        if not faa.is_file():
            raise ValueError(
                f'Missing proteome: {faa}'
            )

        match = ACC_RE.search(str(faa))

        if not match:
            raise ValueError(
                f'Could not extract GCF/GCA accession from: {faa}'
            )

        assembly = match.group(1)

        if assembly in seen_assemblies:
            raise ValueError(
                f'Duplicate genome ID in FAA list: {assembly}'
            )

        seen_assemblies.add(assembly)

        st = faa.stat()

        cohort.append({
            'assembly_accession': assembly,
            'culsans_positive': 'NA'
        })

        inputs.append([
            assembly,
            str(faa),
            st.st_size,
            st.st_mtime_ns
        ])

    if not cohort:
        raise ValueError(
            'FAA list contained no usable proteomes'
        )

    return cohort, inputs


def summarize(out, cohort, info, batches):

    counts = defaultdict(Counter)
    unique = defaultdict(Counter)

    ambiguous = 0
    total_hits = 0

    with (out / 'surface_ko_evidence.tsv').open('w') as f:

        w = csv.writer(f, delimiter='\t')

        w.writerow([
            'assembly',
            'protein',
            'KO',
            'threshold',
            'score',
            'evalue',
            'definition',
            'pathway_membership',
            'multiple_passing_KOs'
        ])

        for batch in batches:

            hits = parse_hits(
                batch / 'kofam_detail.tsv'
            )

            groups = defaultdict(set)

            for r in hits:
                groups[tuple(r[:2])].add(r[2])

            ambiguous += sum(
                len(v) > 1
                for v in groups.values()
            )

            for r in hits:

                amb = (
                    len(groups[tuple(r[:2])]) > 1
                )

                w.writerow(
                    r + [
                        ';'.join(
                            info['mapping'][r[2]]
                        ),
                        'yes' if amb else 'no'
                    ]
                )

                counts[r[0]][r[2]] += 1

                if not amb:
                    unique[r[0]][r[2]] += 1

                total_hits += 1

    outputs = [
        (
            'ko_copy_numbers_all_hits.tsv',
            counts,
            False
        ),
        (
            'ko_copy_numbers_unambiguous.tsv',
            unique,
            False
        ),
        (
            'ko_presence_unambiguous.tsv',
            unique,
            True
        )
    ]

    for name, data, presence in outputs:

        with (out / name).open('w') as f:

            w = csv.writer(f, delimiter='\t')

            w.writerow(
                [
                    'assembly_accession',
                    'culsans_positive'
                ] + info['usable']
            )

            for r in cohort:

                vals = [
                    data[r['assembly_accession']][k]
                    for k in info['usable']
                ]

                if presence:
                    vals = [
                        int(v > 0)
                        for v in vals
                    ]

                w.writerow([
                    r['assembly_accession'],
                    r['culsans_positive']
                ] + vals)

    save(
        out / 'summary.json',
        {
            'genomes': len(cohort),
            'tested_KOs': len(info['usable']),
            'passing_protein_KO_hits': total_hits,
            'proteins_with_multiple_passing_KOs':
                ambiguous,
            'interpretation':
                'Zero means no qualifying detection, not '
                'proven biological absence. Unambiguous '
                'means among the targeted profiles only. '
                'Pathway membership does not prove '
                'surface-structure-specific function, pathway completeness, '
                'or serotype. No association statistics or '
                'phylogenetic correction performed here.'
        }
    )


def main():

    p = argparse.ArgumentParser(
        description=__doc__
    )

    p.add_argument(
        '--root',
        type=Path,
        default=Path(__file__).resolve().parent
    )

    p.add_argument(
        '--faa-list',
        type=Path,
        required=True,
        help='Text file containing one .faa path per genome'
    )

    p.add_argument(
        '--tool-bin',
        type=Path,
        default=Path(
            '/Users/giusym/miniconda3/envs/gtotree/bin'
        )
    )

    p.add_argument(
        '--cpu',
        type=int,
        default=4
    )

    p.add_argument(
        '--batch-size',
        type=int,
        default=20
    )

    p.add_argument(
        '--limit',
        type=int,
        default=0,
        help='Pilot only: number of genomes'
    )

    p.add_argument(
        '--outdir',
        type=Path
    )

    a = p.parse_args()

    if (
        a.cpu < 1
        or a.batch_size < 1
        or a.limit < 0
    ):
        p.error('Invalid numeric argument')

    root = a.root.resolve()

    db = root / 'kofam_database'

    out = (
        a.outdir
        or root / 'results_surface_kofam'
    ).resolve()

    out.mkdir(
        parents=True,
        exist_ok=True
    )

    progress = out / 'progress.json'

    # --------------------------------------------------------
    # Lock output directory
    # --------------------------------------------------------

    import fcntl

    lock = (out / 'run.lock').open('w')

    fcntl.flock(
        lock,
        fcntl.LOCK_EX | fcntl.LOCK_NB
    )

    try:

        # ----------------------------------------------------
        # Prepare KOfam database
        # ----------------------------------------------------

        save(
            progress,
            {
                'stage': 'preparing_database',
                'pid': os.getpid()
            }
        )

        info = prepare(db)

        # ----------------------------------------------------
        # NEW INPUT HANDLING
        # ----------------------------------------------------

        cohort, inputs = read_faa_list(
            a.faa_list
        )

        if a.limit:
            cohort = cohort[:a.limit]
            inputs = inputs[:a.limit]

        print(
            f'Proteomes in cohort: {len(cohort)}',
            flush=True
        )

        print(
            f'Target KOs with calibrated thresholds: '
            f'{len(info["usable"])}',
            flush=True
        )

        # ----------------------------------------------------
        # Manifest
        # ----------------------------------------------------

        manifest = {
            'cohort': cohort,
            'inputs': inputs,
            'database': info,
            'batch_size': a.batch_size,
            'tool_bin':
                str(a.tool_bin.resolve()),
            'script_sha256':
                sha(Path(__file__))
        }

        manifest_path = (
            out / 'manifest.json'
        )

        if manifest_path.exists():

            old_manifest = json.loads(
                manifest_path.read_text()
            )

            if old_manifest != manifest:

                raise ValueError(
                    'Run inputs/settings changed; '
                    'use a fresh output directory'
                )

        save(
            manifest_path,
            manifest
        )

        env = dict(
            os.environ,
            PATH=
                str(a.tool_bin.resolve())
                + os.pathsep
                + os.environ.get('PATH', '')
        )

        # ----------------------------------------------------
        # Run batches
        # ----------------------------------------------------

        batches = []

        total = (
            len(cohort)
            + a.batch_size
            - 1
        ) // a.batch_size

        for i in range(
            0,
            len(cohort),
            a.batch_size
        ):

            batch = (
                out
                / 'batches'
                / f'{i // a.batch_size + 1:04d}'
            )

            batch.mkdir(
                parents=True,
                exist_ok=True
            )

            batches.append(batch)

            done = batch / 'done.json'

            if done.exists():
                continue

            save(
                progress,
                {
                    'stage':
                        'running_kofamscan',
                    'pid':
                        os.getpid(),
                    'batch':
                        len(batches),
                    'total_batches':
                        total,
                    'completed_genomes':
                        i,
                    'total_genomes':
                        len(cohort)
                }
            )

            print(
                f'Batch {len(batches)}/{total}: '
                f'genomes {i + 1}-'
                f'{min(i + a.batch_size, len(cohort))}',
                flush=True
            )

            query = batch / 'proteins.faa'

            nseq = 0
            seen = set()

            with query.open('w') as dest:

                for (
                    assembly,
                    path,
                    _,
                    _
                ) in inputs[
                    i:i + a.batch_size
                ]:

                    per_genome = 0

                    with open(path) as src:

                        for line in src:

                            if line.startswith('>'):

                                protein = (
                                    line[1:]
                                    .split()[0]
                                )

                                name = (
                                    assembly
                                    + '|'
                                    + protein
                                )

                                if name in seen:
                                    raise ValueError(
                                        'Duplicate protein ID: '
                                        + name
                                    )

                                seen.add(name)

                                dest.write(
                                    '>'
                                    + name
                                    + '\n'
                                )

                                nseq += 1
                                per_genome += 1

                            else:
                                dest.write(line)

                    if not per_genome:
                        raise ValueError(
                            'Empty proteome: '
                            + path
                        )

            detail = (
                batch
                / 'kofam_detail.tsv'
            )

            cmd = [
                str(
                    a.tool_bin
                    / 'exec_annotation'
                ),
                '-o',
                str(detail),
                '--profile',
                str(db / 'surface.hal'),
                '--ko-list',
                str(db / 'ko_list'),
                '--cpu',
                str(a.cpu),
                '--tmp-dir',
                str(batch / 'tmp'),
                '--format',
                'detail-tsv',
                str(query)
            ]

            with (
                batch / 'run.log'
            ).open('w') as log:

                subprocess.run(
                    cmd,
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True
                )

            hits = parse_hits(detail)

            save(
                done,
                {
                    'genomes':
                        len(
                            inputs[
                                i:i + a.batch_size
                            ]
                        ),
                    'proteins':
                        nseq,
                    'passing_hits':
                        len(hits)
                }
            )

            # Delete only temporary merged query.
            # Original FAA files remain untouched.
            query.unlink()

        # ----------------------------------------------------
        # Summarise
        # ----------------------------------------------------

        summarize(
            out,
            cohort,
            info,
            batches
        )

        save(
            progress,
            {
                'stage': 'complete',
                'genomes': len(cohort),
                'batches': len(batches)
            }
        )

        print(
            'Completed: ' + str(out),
            flush=True
        )

    except Exception as e:

        save(
            progress,
            {
                'stage': 'failed',
                'error': str(e),
                'pid': os.getpid()
            }
        )

        raise


if __name__ == '__main__':
    main()
