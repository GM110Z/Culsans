#!/usr/bin/env python3
"""Profile defence and anti-defence systems around representative protein hits.

FlankGuard maps protein accessions through NCBI IPG, selects one deterministic
RefSeq nucleotide representative per protein, extracts a genomic neighbourhood,
runs DefenseFinder (including AntiDefenseFinder) and PADLOC, and combines and
visualises the results.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import html
from http.client import IncompleteRead
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from typing import Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
UNKNOWN = "Unknown"
UNKNOWN_VALUES = {
    "", "missing", "unknown", "not applicable", "not_applicable", "n/a",
    "na", "none", "not provided", "not collected", "not determined", "nd",
    "null", "nan", "-", "--",
}
ACCESSION_RE = re.compile(r"(?:[A-Z]{1,4}_)?[A-Z]{0,3}\d{5,}(?:\.\d+)?", re.I)


def clean(value: object) -> str:
    if value is None:
        return UNKNOWN
    text = str(value).strip()
    return UNKNOWN if text.lower().rstrip(".") in UNKNOWN_VALUES else text


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1].lower()


def normalise_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def attr(element: ET.Element, *names: str) -> str:
    attributes = {normalise_key(k): v for k, v in element.attrib.items()}
    for name in names:
        value = clean(attributes.get(normalise_key(name)))
        if value != UNKNOWN:
            return value
    return UNKNOWN


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", value)


def versionless(value: str) -> str:
    return clean(value).split(".", 1)[0]


class NcbiClient:
    def __init__(self, email: str, api_key: str | None, delay: float, retries: int = 8):
        self.email = email
        self.api_key = api_key
        self.delay = delay
        self.retries = retries
        self.last_request = 0.0

    def get(self, endpoint: str, params: dict[str, object]) -> str:
        query: dict[str, object] = {
            **params, "email": self.email, "tool": "flankguard",
        }
        if self.api_key:
            query["api_key"] = self.api_key
        url = BASE + endpoint + "?" + urlencode(query)
        for attempt in range(self.retries):
            pause = self.delay - (time.monotonic() - self.last_request)
            if pause > 0:
                time.sleep(pause)
            request = Request(
                url,
                headers={"User-Agent": f"FlankGuard/1.0 ({self.email})"},
            )
            try:
                with urlopen(request, timeout=120) as response:
                    result = response.read().decode("utf-8")
                self.last_request = time.monotonic()
                return result
            except HTTPError as exc:
                self.last_request = time.monotonic()
                if exc.code not in {429, 500, 502, 503, 504}:
                    raise RuntimeError(f"NCBI HTTP {exc.code} for {endpoint}") from exc
                reason = f"HTTP {exc.code}"
            except (URLError, TimeoutError, IncompleteRead) as exc:
                self.last_request = time.monotonic()
                reason = str(exc)
            wait = max(1.0, self.delay) * (2**attempt)
            print(f"  NCBI request failed ({reason}); retrying in {wait:.1f}s", file=sys.stderr)
            time.sleep(wait)
        raise RuntimeError(f"NCBI request failed after {self.retries} attempts: {endpoint}")


def cache_fetch(
    client: NcbiClient,
    path: Path,
    endpoint: str,
    params: dict[str, object],
    offline: bool,
) -> str:
    if path.exists() and path.stat().st_size:
        return path.read_text(encoding="utf-8")
    if offline:
        raise FileNotFoundError(f"Required cache file is absent in offline mode: {path}")
    text = client.get(endpoint, params)
    atomic_write(path, text)
    return text


def read_proteins(path: Path, protein_column: str | None) -> tuple[list[str], int, str]:
    sample = path.read_text(encoding="utf-8-sig", errors="replace")[:65536]
    if not sample.strip():
        raise ValueError(f"Input is empty: {path}")
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters=",\t;|").delimiter
    except csv.Error:
        delimiter = "\t" if "\t" in sample.partition("\n")[0] else ","
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
        raw_rows = list(csv.reader(handle, delimiter=delimiter))
    if not raw_rows or len(raw_rows[0]) < 2:
        raise ValueError("Input table must contain at least two columns")
    first = [cell.strip() for cell in raw_rows[0]]
    normalised = {normalise_key(name): index for index, name in enumerate(first)}
    recognised = [
        "proteinaccession", "sseqid", "subjectaccession", "subjectid", "accession",
    ]
    if protein_column:
        wanted = normalise_key(protein_column)
        if wanted not in normalised:
            raise ValueError(f"Protein column not found: {protein_column}")
        index = normalised[wanted]
        has_header = True
        column_name = first[index]
    else:
        matches = [name for name in recognised if name in normalised]
        if matches:
            index = normalised[matches[0]]
            has_header = True
            column_name = first[index]
        else:
            index = 1
            has_header = False
            column_name = "column_2"
    data = raw_rows[1:] if has_header else raw_rows
    proteins: list[str] = []
    seen: set[str] = set()
    for row in data:
        if len(row) <= index:
            continue
        matches = ACCESSION_RE.findall(row[index])
        accession = (matches[-1] if matches else row[index].split()[0]).upper() if row[index].strip() else ""
        if accession and accession not in seen:
            proteins.append(accession)
            seen.add(accession)
    if not proteins:
        raise ValueError(f"No protein accessions found in {column_name}")
    return proteins, len(data), column_name


def parse_ipg_xml(xml_text: str, requested: str) -> list[dict[str, str]]:
    root = ET.fromstring(xml_text)
    mappings: list[dict[str, str]] = []

    def walk(element: ET.Element, protein: dict[str, str] | None = None) -> None:
        tag = local_name(element.tag)
        if tag in {"protein", "proteinrecord"}:
            protein = {
                "ipg_source": attr(element, "source", "database", "db"),
                "ipg_protein_accession": attr(element, "accver", "accession", "protein_acc"),
                "protein_name": attr(element, "name", "product", "protein_name"),
            }
        if tag in {"cds", "cdsrecord", "nucleotide"}:
            source = (protein or {}).get("ipg_source", UNKNOWN)
            if source == UNKNOWN:
                source = attr(element, "source", "database", "db")
            nucleotide = attr(
                element, "accver", "nucleotide_acc", "nucleotide_accession",
                "nucaccver", "nuc_accession", "accession",
            )
            if source.lower() == "refseq" and nucleotide != UNKNOWN:
                mappings.append({
                    "protein_accession": requested,
                    "ipg_protein_accession": (protein or {}).get("ipg_protein_accession", requested),
                    "protein_name": (protein or {}).get("protein_name", UNKNOWN),
                    "nuccore_accession": nucleotide,
                    "cds_start": attr(element, "start", "from"),
                    "cds_stop": attr(element, "stop", "to", "end"),
                    "strand": attr(element, "strand"),
                    "organism": attr(element, "organism", "org"),
                    "strain": attr(element, "strain"),
                    "assembly": attr(element, "assembly", "assembly_acc", "assembly_accession"),
                })
        for child in element:
            walk(child, protein)

    walk(root)
    unique: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in mappings:
        key = (row["nuccore_accession"], row["cds_start"], row["cds_stop"])
        unique[key] = row
    return list(unique.values())


def accession_rank(accession: str) -> int:
    upper = accession.upper()
    if upper.startswith("NC_"):
        return 0
    if upper.startswith("NZ_CP"):
        return 1
    if upper.startswith(("NZ_AP", "NZ_BA")):
        return 2
    if upper.startswith("NZ_"):
        return 3
    return 4


def valid_coordinate(value: str) -> bool:
    return bool(re.fullmatch(r"\d+", clean(value)))


def choose_representative(rows: list[dict[str, str]]) -> dict[str, str]:
    if not rows:
        raise ValueError("No RefSeq IPG mappings were available")

    def ranking(row: dict[str, str]) -> tuple[object, ...]:
        valid = valid_coordinate(row["cds_start"]) and valid_coordinate(row["cds_stop"])
        return (
            0 if valid else 1,
            accession_rank(row["nuccore_accession"]),
            0 if clean(row["assembly"]) != UNKNOWN else 1,
            0 if clean(row["organism"]) != UNKNOWN else 1,
            0 if clean(row["strain"]) != UNKNOWN else 1,
            row["nuccore_accession"],
        )

    selected = min(rows, key=ranking).copy()
    selected["ipg_refseq_candidates"] = str(len(rows))
    selected["representative_rank"] = str(accession_rank(selected["nuccore_accession"]))
    selected["representative_rule"] = (
        "valid CDS coordinates; prefer NC_, then NZ_CP, NZ_AP/NZ_BA, other NZ_; "
        "prefer assembly, organism and strain metadata; lexical tie-break"
    )
    return selected


def parse_nuccore_summary(xml_text: str) -> dict[str, str]:
    root = ET.fromstring(xml_text)
    values: dict[str, str] = {}
    for element in root.iter():
        if local_name(element.tag) != "item":
            continue
        name = clean(element.attrib.get("Name") or element.attrib.get("name"))
        value = clean(element.text)
        if name != UNKNOWN and value != UNKNOWN:
            values[normalise_key(name)] = value
    return {
        "replicon_length": values.get("length", UNKNOWN),
        "replicon_topology": values.get("topology", UNKNOWN),
        "replicon_title": values.get("title", UNKNOWN),
    }


def requested_intervals(
    cds_start: int,
    cds_stop: int,
    flank: int,
    length: int | None,
    topology: str,
) -> tuple[list[tuple[int, int]], str]:
    desired_start = cds_start - flank
    desired_stop = cds_stop + flank
    if length is None or length < 1:
        return [(max(1, desired_start), desired_stop)], "length unknown; ends may be clamped by NCBI"
    if topology.lower() == "circular" and desired_stop - desired_start + 1 < length:
        start = ((desired_start - 1) % length) + 1
        stop = ((desired_stop - 1) % length) + 1
        if desired_start < 1 or desired_stop > length:
            return [(start, length), (1, stop)], "complete circular wrap"
        return [(desired_start, desired_stop)], "complete"
    start = max(1, desired_start)
    stop = min(length, desired_stop)
    status = "complete" if start == desired_start and stop == desired_stop else "truncated at linear replicon edge"
    return [(start, stop)], status


def fetch_context(
    client: NcbiClient,
    representative: dict[str, str],
    flank: int,
    cache: Path,
    offline: bool,
) -> tuple[str, str, list[dict[str, object]], dict[str, str]]:
    accession = representative["nuccore_accession"]
    summary_xml = cache_fetch(
        client, cache / "nuccore_summary" / f"{safe_name(accession)}.xml",
        "esummary.fcgi", {"db": "nuccore", "id": accession, "retmode": "xml"}, offline,
    )
    summary = parse_nuccore_summary(summary_xml)
    length = int(summary["replicon_length"]) if valid_coordinate(summary["replicon_length"]) else None
    cds_start = min(int(representative["cds_start"]), int(representative["cds_stop"]))
    cds_stop = max(int(representative["cds_start"]), int(representative["cds_stop"]))
    intervals, edge_status = requested_intervals(
        cds_start, cds_stop, flank, length, summary["replicon_topology"],
    )
    all_sequence: list[str] = []
    all_features: list[dict[str, object]] = []
    piece_metadata: list[dict[str, str]] = []
    returned_accession = accession
    offset = 0
    for start, stop in intervals:
        region_key = f"{safe_name(accession)}_{start}_{stop}.xml"
        region_xml = cache_fetch(client, cache / "regions" / region_key, "efetch.fcgi", {
            "db": "nuccore", "id": accession, "rettype": "gbwithparts", "retmode": "xml",
            "seq_start": start, "seq_stop": stop,
        }, offline)
        piece_accession, sequence, features, metadata = parse_region_xml(
            region_xml, representative, start, stop,
        )
        if piece_accession != UNKNOWN:
            returned_accession = piece_accession
        for feature in features:
            feature["local_start"] = int(feature["local_start"]) + offset
            feature["local_stop"] = int(feature["local_stop"]) + offset
            all_features.append(feature)
        all_sequence.append(sequence)
        piece_metadata.append(metadata)
        offset += len(sequence)
    for index, feature in enumerate(all_features, 1):
        feature["feature_id"] = f"FG_{index:05d}"
    metadata = {
        "returned_accession": returned_accession,
        "definition": piece_metadata[0]["definition"],
        "organism_from_record": piece_metadata[0]["organism_from_record"],
        "request_start": str(intervals[0][0]),
        "requested_stop": str(intervals[-1][1]),
        "request_intervals": ";".join(f"{start}-{stop}" for start, stop in intervals),
        "returned_length_bp": str(sum(len(sequence) for sequence in all_sequence)),
        "returned_stop": str(intervals[-1][1]),
        "coordinate_mode": ";".join(dict.fromkeys(item["coordinate_mode"] for item in piece_metadata)),
        "target_cds_found": "yes" if any(item["target_cds_found"] == "yes" for item in piece_metadata) else "no",
        "cds_count": str(len(all_features)),
        "replicon_length": summary["replicon_length"],
        "replicon_topology": summary["replicon_topology"],
        "edge_status": edge_status,
    }
    return returned_accession, "".join(all_sequence), all_features, metadata


def element_text(parent: ET.Element, suffix: str) -> str:
    for element in parent.iter():
        if local_name(element.tag).endswith(suffix.lower()):
            value = clean(element.text)
            if value != UNKNOWN:
                return value
    return UNKNOWN


def feature_qualifiers(feature: ET.Element) -> dict[str, list[str]]:
    result: dict[str, list[str]] = defaultdict(list)
    for qualifier in feature.iter():
        if local_name(qualifier.tag) not in {"gbqualifier", "insdqualifier"}:
            continue
        name = UNKNOWN
        value = UNKNOWN
        for child in qualifier:
            tag = local_name(child.tag)
            if tag.endswith("_name"):
                name = clean(child.text)
            elif tag.endswith("_value"):
                value = clean(child.text)
        if name != UNKNOWN and value != UNKNOWN:
            result[name.lower()].append(value)
    return result


def first_value(mapping: dict[str, list[str]], *names: str) -> str:
    for name in names:
        values = mapping.get(name.lower(), [])
        if values:
            return clean(values[0])
    return UNKNOWN


def parse_region_xml(
    xml_text: str,
    representative: dict[str, str],
    request_start: int,
    requested_stop: int,
) -> tuple[str, str, list[dict[str, object]], dict[str, str]]:
    root = ET.fromstring(xml_text)
    record = next(
        (element for element in root.iter() if local_name(element.tag) in {"gbseq", "insdseq"}),
        None,
    )
    if record is None:
        raise ValueError("No GenBank sequence record in NCBI response")
    sequence = element_text(record, "_sequence").replace(" ", "").replace("\n", "").upper()
    if sequence == UNKNOWN or not sequence:
        raise ValueError("NCBI region response contained no nucleotide sequence")
    accession = element_text(record, "_accession-version")
    definition = element_text(record, "_definition")
    organism = element_text(record, "_organism")
    features: list[dict[str, object]] = []
    raw_intervals: list[tuple[int, int, str, ET.Element]] = []
    for feature in record.iter():
        if local_name(feature.tag) not in {"gbfeature", "insdfeature"}:
            continue
        key = element_text(feature, "_key")
        if key.lower() != "cds":
            continue
        intervals = []
        strand = "+"
        for interval in feature.iter():
            if local_name(interval.tag) not in {"gbinterval", "insdinterval"}:
                continue
            start_text = element_text(interval, "_from")
            stop_text = element_text(interval, "_to")
            point_text = element_text(interval, "_point")
            if valid_coordinate(start_text) and valid_coordinate(stop_text):
                start_zero, stop_zero = int(start_text), int(stop_text)
            elif valid_coordinate(point_text):
                start_zero = stop_zero = int(point_text)
            else:
                continue
            is_comp = element_text(interval, "_iscomp")
            if is_comp == UNKNOWN:
                for child in interval:
                    if local_name(child.tag).endswith("_iscomp"):
                        is_comp = attr(child, "value")
                        break
            if is_comp.lower() in {"true", "1", "yes"}:
                strand = "-"
            intervals.append((min(start_zero, stop_zero), max(start_zero, stop_zero)))
        if intervals:
            raw_intervals.append((min(x[0] for x in intervals), max(x[1] for x in intervals), strand, feature))
    if not raw_intervals:
        raise ValueError("No translated CDS features were returned for the selected region")

    # For an EFetch request with seq_start/seq_stop, GBInterval coordinates are
    # 1-based and relative to the returned slice (the feature location strings
    # use the same coordinates). Convert them back to genomic coordinates for
    # the manifest while retaining slice coordinates for FAA/GFF generation.
    target_accessions = {
        versionless(representative["protein_accession"]),
        versionless(representative["ipg_protein_accession"]),
    }
    target_start = min(int(representative["cds_start"]), int(representative["cds_stop"]))
    target_stop = max(int(representative["cds_start"]), int(representative["cds_stop"]))
    target_found = False
    for raw_start, raw_stop, strand, feature in raw_intervals:
        local_start, local_stop = raw_start, raw_stop
        genomic_start = request_start + raw_start - 1
        genomic_stop = request_start + raw_stop - 1
        local_start = max(1, local_start)
        local_stop = min(len(sequence), local_stop)
        if local_stop < local_start:
            continue
        qualifiers = feature_qualifiers(feature)
        translation = first_value(qualifiers, "translation").replace(" ", "").replace("\n", "")
        if translation == UNKNOWN or not translation:
            continue
        protein_id = first_value(qualifiers, "protein_id", "locus_tag", "gene")
        locus_tag = first_value(qualifiers, "locus_tag")
        product = first_value(qualifiers, "product")
        by_accession = versionless(protein_id) in target_accessions
        by_coordinates = genomic_start <= target_stop and genomic_stop >= target_start
        is_target = by_accession or by_coordinates
        target_found = target_found or is_target
        features.append({
            "local_start": local_start,
            "local_stop": local_stop,
            "genomic_start": genomic_start,
            "genomic_stop": genomic_stop,
            "strand": strand,
            "original_protein_id": protein_id,
            "locus_tag": locus_tag,
            "product": product,
            "translation": translation,
            "is_target": "yes" if is_target else "no",
        })
    features.sort(key=lambda row: (int(row["local_start"]), int(row["local_stop"])))
    if not features:
        raise ValueError("No translated CDS features could be extracted from the region")
    for index, feature in enumerate(features, 1):
        feature["feature_id"] = f"FG_{index:05d}"
    metadata = {
        "returned_accession": accession,
        "definition": definition,
        "organism_from_record": organism,
        "request_start": str(request_start),
        "requested_stop": str(requested_stop),
        "returned_length_bp": str(len(sequence)),
        "returned_stop": str(request_start + len(sequence) - 1),
        "coordinate_mode": "slice-relative-1-based",
        "target_cds_found": "yes" if target_found else "no",
        "cds_count": str(len(features)),
    }
    return accession, sequence, features, metadata


def wrap_fasta(sequence: str, width: int = 70) -> str:
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def write_context_files(
    directory: Path,
    representative: dict[str, str],
    accession: str,
    sequence: str,
    features: list[dict[str, object]],
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    contig = safe_name(f"{representative['protein_accession']}__{accession}")
    atomic_write(directory / "neighborhood.fna", f">{contig}\n{wrap_fasta(sequence)}\n")
    faa_parts = []
    gff_parts = ["##gff-version 3"]
    gene_rows = []
    for feature in features:
        feature_id = str(feature["feature_id"])
        product = str(feature["product"])
        faa_parts.append(f">{feature_id} {product}\n{wrap_fasta(str(feature['translation']))}\n")
        attributes = (
            f"ID={feature_id};protein_id={feature_id};"
            f"original_protein_id={safe_name(str(feature['original_protein_id']))};"
            f"product={safe_name(product)}"
        )
        gff_parts.append("\t".join([
            contig, "FlankGuard", "CDS", str(feature["local_start"]),
            str(feature["local_stop"]), ".", str(feature["strand"]), "0", attributes,
        ]))
        gene_rows.append({key: feature.get(key, "") for key in [
            "feature_id", "local_start", "local_stop", "genomic_start", "genomic_stop",
            "strand", "original_protein_id", "locus_tag", "product", "is_target",
        ]})
    atomic_write(directory / "neighborhood.faa", "".join(faa_parts))
    atomic_write(directory / "neighborhood.gff", "\n".join(gff_parts) + "\n")
    write_tsv(directory / "genes.tsv", gene_rows, [
        "feature_id", "local_start", "local_stop", "genomic_start", "genomic_stop",
        "strand", "original_protein_id", "locus_tag", "product", "is_target",
    ])


def conda_prefix(conda: str, environment: str) -> list[str]:
    # Some desktop/IDE launchers prepend system paths ahead of the environment
    # even inside `conda run`. Force the selected environment's bin directory
    # to the front so tools such as PADLOC use their matching R/Rscript.
    return [
        conda, "run", "--no-capture-output", "-n", environment,
        "/bin/bash", "-c", 'export PATH="$CONDA_PREFIX/bin:$PATH"; exec "$@"',
        "flankguard-environment",
    ]


def run_logged(command: list[str], log: Path, dry_run: bool) -> None:
    printable = shlex.join(command)
    print(f"  $ {printable}", flush=True)
    if dry_run:
        atomic_write(log, "DRY RUN\n" + printable + "\n")
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        handle.write("COMMAND\n" + printable + "\n\nOUTPUT\n")
        handle.flush()
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if completed.returncode:
        raise RuntimeError(f"Command failed with exit code {completed.returncode}; see {log}")


def check_tool(command: list[str], label: str, reject_missing_warning: bool = False) -> str:
    completed = subprocess.run(command, capture_output=True, text=True)
    output = (completed.stdout or "") + (completed.stderr or "")
    if completed.returncode:
        raise RuntimeError(
            f"{label} preflight check failed (exit {completed.returncode}):\n{output.strip()}"
        )
    if reject_missing_warning and re.search(r"warning:.*not installed", output, re.I):
        missing = "; ".join(
            line.strip() for line in output.splitlines()
            if re.search(r"warning:.*not installed", line, re.I)
        )
        raise RuntimeError(f"{label} dependency check failed: {missing}")
    return output


def locate_files(directory: Path, patterns: tuple[str, ...]) -> list[Path]:
    found: list[Path] = []
    for pattern in patterns:
        found.extend(directory.rglob(pattern))
    return sorted(set(found))


def read_delimited(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def value_from(row: dict[str, str], *candidates: str) -> str:
    indexed = {normalise_key(key): clean(value) for key, value in row.items()}
    for candidate in candidates:
        value = indexed.get(normalise_key(candidate), UNKNOWN)
        if value != UNKNOWN:
            return value
    return UNKNOWN


def parse_defensefinder(directory: Path) -> list[dict[str, str]]:
    results = []
    files = locate_files(directory, ("*_defense_finder_systems.tsv", "*systems.tsv"))
    for path in files:
        for row in read_delimited(path, "\t"):
            combined = " ".join(str(value) for value in row.values()).lower()
            category = "AntiDefenseFinder" if "antidefense" in combined or "anti-defense" in combined else "DefenseFinder"
            system = value_from(row, "subtype", "type", "system", "system_name", "model_fqn", "model")
            system_id = value_from(row, "sys_id", "system_id", "id")
            results.append({
                "tool": category,
                "system": system,
                "system_id": system_id,
                "result_file": str(path),
                "raw_result": json.dumps(row, sort_keys=True),
            })
    return results


def parse_padloc(directory: Path) -> list[dict[str, str]]:
    results = []
    files = locate_files(directory, ("*_padloc.csv",))
    for path in files:
        rows = read_delimited(path, ",")
        grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            system = value_from(row, "system", "system_name")
            number = value_from(row, "system.number", "system_number", "system.id")
            grouped[(system, number)].append(row)
        for (system, number), members in grouped.items():
            results.append({
                "tool": "PADLOC",
                "system": system,
                "system_id": number,
                "result_file": str(path),
                "raw_result": json.dumps(members, sort_keys=True),
            })
    return results


def hex_to_rgb(colour: str) -> tuple[int, int, int]:
    value = colour.lstrip("#")
    if not re.fullmatch(r"[0-9A-Fa-f]{6}", value):
        raise ValueError(f"Colour must be a six-digit hexadecimal value: {colour}")
    return tuple(int(value[index:index + 2], 16) for index in (0, 2, 4))  # type: ignore[return-value]


def tint(colour: str, fraction: float) -> str:
    red, green, blue = hex_to_rgb(colour)
    mixed = tuple(round(channel + (255 - channel) * fraction) for channel in (red, green, blue))
    return "#" + "".join(f"{channel:02X}" for channel in mixed)


def make_heatmap(
    summary_rows: list[dict[str, object]],
    presence: dict[str, set[str]],
    frequencies: Counter[str],
    colour: str,
    top_systems: int,
    max_rows: int,
) -> str:
    systems = [name for name, _ in frequencies.most_common(top_systems)]
    selected = set(systems)
    eligible = [
        row for row in summary_rows
        if presence.get(str(row["protein_accession"]), set()) & selected
    ]
    ranked = sorted(
        eligible,
        key=lambda row: (
            -len(presence.get(str(row["protein_accession"]), set()) & selected),
            tuple(-int(system in presence.get(str(row["protein_accession"]), set())) for system in systems),
            str(row["organism"]),
            str(row["protein_accession"]),
        ),
    )[:max_rows]
    cell = 20
    left = 430
    top = 285
    width = max(900, left + cell * max(1, len(systems)) + 120)
    height = top + cell * max(1, len(ranked)) + 95
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>.title{{font:700 25px Arial,sans-serif;fill:#111}}.note{{font:13px Arial,sans-serif;fill:#444}}.rowlabel{{font:11px Arial,sans-serif;fill:#111}}.collabel{{font:12px Arial,sans-serif;fill:#111}}</style>
<rect width="100%" height="100%" fill="#FFFFFF"/>
<text x="35" y="42" class="title">Defence-system presence around Culsans</text>
<text x="35" y="68" class="note">Rows are representative RefSeq neighbourhoods; columns are the most frequent tool-specific system calls.</text>
<text x="35" y="89" class="note">Blue = detected; white = not detected. Rows without any displayed top-system call are omitted.</text>''']
    if not systems or not ranked:
        parts.append('<text x="35" y="115" class="note">No systems were detected.</text></svg>')
        return "\n".join(parts)
    for column, system in enumerate(systems):
        x = left + column * cell + 15
        display = system.replace("DefenseFinder: ", "DF: ").replace("AntiDefenseFinder: ", "ADF: ").replace("_", " ")
        parts.append(f'<text x="{x}" y="{top - 12}" transform="rotate(-58 {x} {top - 12})" class="collabel">{html.escape(display)} (n={frequencies[system]})</text>')
    for row_number, row in enumerate(ranked):
        y = top + row_number * cell
        protein = str(row["protein_accession"])
        organism = str(row.get("organism", "")).strip() or "Unknown organism"
        parts.extend([
            f'<text x="{left - 15}" y="{y + 14}" text-anchor="end" class="rowlabel"><tspan font-style="italic">{html.escape(organism)}</tspan><tspan font-style="normal"> | {html.escape(protein)}</tspan></text>',
            f'<line x1="{left - 12}" y1="{y + 10}" x2="{left}" y2="{y + 10}" stroke="#111" stroke-width="1"/>',
        ])
        for column, system in enumerate(systems):
            x = left + column * cell
            fill = colour if system in presence.get(protein, set()) else "#FFFFFF"
            parts.append(f'<rect x="{x}" y="{y}" width="20" height="20" fill="{fill}" stroke="#111" stroke-width="1"/>')
    parts.append(f'<text x="35" y="{height - 48}" class="note">Showing {len(ranked)} of {len(eligible)} neighbourhoods containing at least one displayed system; {len(systems)} system columns.</text>')
    parts.append(f'<rect x="35" y="{height - 31}" width="14" height="14" fill="{colour}" stroke="#111"/><text x="57" y="{height - 19}" class="note">Detected</text>')
    parts.append(f'<rect x="145" y="{height - 31}" width="14" height="14" fill="#FFFFFF" stroke="#111"/><text x="167" y="{height - 19}" class="note">Not detected</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def make_frequency_svg(rows: list[dict[str, object]], colour: str, top: int, tool_specific: bool = True) -> str:
    shown = rows[:top]
    width = 1050
    row_height = 35
    plot_left = 385
    plot_width = 500
    plot_top = 120
    height = plot_top + row_height * max(1, len(shown)) + 105
    percent_field = "percent_of_neighborhoods" if tool_specific else "percent_of_analysed_neighborhoods"
    maximum = max([float(row[percent_field]) for row in shown] or [1.0])
    axis_max = max(5.0, float(math.ceil(maximum / 5.0) * 5))
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>.title{{font:700 25px Arial,sans-serif;fill:#111}}.note{{font:13px Arial,sans-serif;fill:#444}}.label{{font:14px Arial,sans-serif;fill:#111}}.value{{font:700 13px Arial,sans-serif;fill:#111}}.tick{{font:12px Arial,sans-serif;fill:#444}}</style>
<rect width="100%" height="100%" fill="#FFFFFF"/>
<text x="40" y="42" class="title">Defence systems associated with Culsans</text>
<text x="40" y="67" class="note">{html.escape('Nonredundant overlapping calls from DefenseFinder and PADLOC are combined by genomic locus.' if not tool_specific else 'Tool-specific calls are shown separately for provenance.')}</text>''']
    if not shown:
        parts.append('<text x="40" y="105" class="note">No systems were detected.</text></svg>')
        return "\n".join(parts)
    for tick in range(0, int(axis_max) + 1, 5):
        x = plot_left + plot_width * tick / axis_max
        parts.append(f'<line x1="{x:.2f}" y1="{plot_top - 12}" x2="{x:.2f}" y2="{plot_top + row_height * len(shown) - 8}" stroke="#D9DEE5" stroke-width="1"/>')
        parts.append(f'<text x="{x:.2f}" y="{plot_top + row_height * len(shown) + 16}" text-anchor="middle" class="tick">{tick}%</text>')
    for index, row in enumerate(shown):
        y = plot_top + index * row_height
        system_field = "system" if tool_specific else "consensus_system"
        display_system = str(row[system_field]).replace("_", " ")
        label = f"{row['tool']}: {display_system}" if tool_specific else display_system
        value = int(row["neighborhood_count"])
        percent = float(row[percent_field])
        point_x = plot_left + plot_width * percent / axis_max
        parts.extend([
            f'<text x="{plot_left - 18}" y="{y + 5}" text-anchor="end" class="label">{html.escape(label)}</text>',
            f'<rect x="{plot_left}" y="{y - 10}" width="{max(2.0, point_x - plot_left):.2f}" height="20" fill="{colour}"/>',
            f'<text x="{point_x + 14:.2f}" y="{y + 5}" class="value">{value} ({percent:.2f}%)</text>',
        ])
    axis_y = plot_top + row_height * len(shown) - 8
    parts.append(f'<line x1="{plot_left}" y1="{axis_y}" x2="{plot_left + plot_width}" y2="{axis_y}" stroke="#111" stroke-width="1.5"/>')
    parts.append(f'<text x="{plot_left + plot_width / 2}" y="{height - 28}" text-anchor="middle" class="label">Culsans neighbourhoods containing system (%)</text>')
    parts.append("</svg>")
    return "\n".join(parts)


def make_frequency_heatmap_svg(rows: list[dict[str, object]], colour: str, top: int) -> str:
    """Draw a compact, one-column system-frequency heatmap."""
    shown = rows[:top]
    row_height = 34
    left = 355
    top_margin = 105
    width = 720
    height = top_margin + row_height * max(1, len(shown)) + 70
    maximum = max([int(row["neighborhood_count"]) for row in shown] or [1])
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>.title{{font:700 24px Arial,sans-serif;fill:#111}}.note{{font:13px Arial,sans-serif;fill:#444}}.label{{font:17px Arial,sans-serif;fill:#111}}.value{{font:14px Arial,sans-serif;fill:#111}}</style>
<rect width="100%" height="100%" fill="#FFFFFF"/>
<text x="24" y="34" class="title">Systems near Culsans</text>
<text x="24" y="58" class="note">Cell intensity represents frequency among analysed RefSeq neighbourhoods (n=483).</text>
<text x="{left + 16}" y="88" text-anchor="middle" class="note">Frequency</text>''']
    if not shown:
        parts.append('<text x="24" y="110" class="note">No systems were detected.</text></svg>')
        return "\n".join(parts)
    for index, row in enumerate(shown):
        y = top_margin + index * row_height
        value = int(row["neighborhood_count"])
        percent = str(row["percent_of_neighborhoods"])
        # Square-root scaling retains visible contrast when the leading system is dominant.
        intensity = 0.12 + 0.88 * (value / maximum) ** 0.5
        fill = tint(colour, 1.0 - intensity)
        label = f"{row['tool']}: {row['system']}"
        parts.extend([
            f'<text x="{left - 14}" y="{y + 23}" text-anchor="end" class="label">{html.escape(label)}</text>',
            f'<line x1="{left - 12}" y1="{y + 17}" x2="{left}" y2="{y + 17}" stroke="#111" stroke-width="2"/>',
            f'<rect x="{left}" y="{y}" width="34" height="34" fill="{fill}" stroke="#111" stroke-width="2"/>',
            f'<text x="{left + 48}" y="{y + 22}" class="value">{value} ({percent}%)</text>',
        ])
    parts.append(f'<text x="24" y="{height - 24}" class="note">Colour endpoint: {html.escape(colour)}. Tool-specific labels are retained.</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def make_a4_taxon_heatmap(
    summary_rows: list[dict[str, object]],
    presence: dict[str, set[str]],
    frequencies: Counter[str],
    colour: str,
    top_systems: int = 15,
    top_genera: int = 12,
) -> str:
    """Create an A4-landscape genus-by-system prevalence heatmap."""
    systems = [name for name, _ in frequencies.most_common(top_systems)]
    genus_rows: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in summary_rows:
        organism = str(row.get("organism", "")).strip()
        genus = organism.split()[0] if organism else "Unknown"
        genus_rows[genus].append(row)
    selected_systems = set(systems)
    informative_genera = [
        genus for genus, members in genus_rows.items()
        if any(presence.get(str(member["protein_accession"]), set()) & selected_systems for member in members)
    ]
    ranked_genera = sorted(informative_genera, key=lambda g: (-len(genus_rows[g]), g))
    displayed = ranked_genera[:top_genera]
    groups: list[tuple[str, list[dict[str, object]]]] = [(g, genus_rows[g]) for g in displayed]
    displayed_set = set(displayed)
    remainder = [row for genus, members in genus_rows.items() if genus not in displayed_set for row in members]
    if remainder:
        groups.append(("Other genera", remainder))

    width, height = 1123, 794  # A4 landscape aspect ratio at screen resolution.
    left, top, cell_w, cell_h = 235, 290, 43, 29
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>.title{{font:700 24px Arial,sans-serif;fill:#111}}.note{{font:13px Arial,sans-serif;fill:#444}}.rowlabel{{font:15px Arial,sans-serif;fill:#111}}.collabel{{font:12px Arial,sans-serif;fill:#111}}.cellvalue{{font:10px Arial,sans-serif}}</style>
<rect width="100%" height="100%" fill="#FFFFFF"/>
<text x="38" y="42" class="title">Defence systems associated with Culsans across bacterial genera</text>
<text x="38" y="69" class="note">Cell values are percentages of representative Culsans neighbourhoods within each genus; n is the number analysed.</text>
<text x="38" y="91" class="note">The 12 best-represented genera with a displayed call are shown separately; remaining genera are pooled.</text>''']
    for column, system in enumerate(systems):
        x = left + column * cell_w + 22
        display = system.replace("DefenseFinder: ", "DF: ").replace("AntiDefenseFinder: ", "ADF: ").replace("_", " ")
        parts.append(f'<text x="{x}" y="{top - 12}" transform="rotate(-58 {x} {top - 12})" class="collabel">{html.escape(display)}</text>')
    for row_index, (genus, members) in enumerate(groups):
        y = top + row_index * cell_h
        label = f"{genus} (n={len(members)})"
        style = '' if genus == "Other genera" else ' font-style="italic"'
        parts.append(f'<text x="{left - 14}" y="{y + 19}" text-anchor="end" class="rowlabel"{style}>{html.escape(label)}</text>')
        proteins = [str(member["protein_accession"]) for member in members]
        for column, system in enumerate(systems):
            percent = 100 * sum(system in presence.get(protein, set()) for protein in proteins) / len(proteins)
            fill = "#FFFFFF" if percent == 0 else tint(colour, 1.0 - percent / 100)
            text_colour = "#FFFFFF" if percent >= 55 else "#111111"
            x = left + column * cell_w
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{fill}" stroke="#111" stroke-width="1"/>')
            if percent > 0:
                parts.append(f'<text x="{x + cell_w / 2}" y="{y + 19}" text-anchor="middle" class="cellvalue" fill="{text_colour}">{percent:.0f}</text>')
    legend_y = top + len(groups) * cell_h + 42
    parts.append(f'<text x="{left}" y="{legend_y}" class="note">Within-genus prevalence (%)</text>')
    for index, value in enumerate((0, 25, 50, 75, 100)):
        x = left + 170 + index * 58
        fill = "#FFFFFF" if value == 0 else tint(colour, 1.0 - value / 100)
        parts.append(f'<rect x="{x}" y="{legend_y - 17}" width="28" height="18" fill="{fill}" stroke="#111"/><text x="{x + 34}" y="{legend_y - 3}" class="note">{value}</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def make_association_breadth_svg(
    summary_rows: list[dict[str, object]],
    consensus_loci: list[dict[str, object]],
    colour: str,
) -> str:
    positive = [row for row in summary_rows if row["defence_environment"] == "yes"]
    total = len(summary_rows)
    associated = len(positive)
    percent = 100 * associated / total if total else 0.0
    genera = {str(row.get("organism", "")).split()[0] for row in positive if str(row.get("organism", "")).strip()}
    locus_counts = Counter(str(row["Culsans_protein_accession"]) for row in consensus_loci)
    multiple = sum(locus_counts[str(row["protein_accession"])] >= 2 for row in positive)
    multiple_percent = 100 * multiple / associated if associated else 0.0
    width, height = 1000, 360
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<style>.title{{font:700 25px Arial,sans-serif;fill:#111}}.big{{font:700 62px Arial,sans-serif;fill:{colour}}}.statement{{font:700 24px Arial,sans-serif;fill:#111}}.metric{{font:700 29px Arial,sans-serif;fill:{colour}}}.label{{font:15px Arial,sans-serif;fill:#333}}.note{{font:13px Arial,sans-serif;fill:#555}}</style>
<rect width="100%" height="100%" fill="#FFFFFF"/>
<text x="45" y="42" class="title">Breadth of Culsans association with defence systems</text>
<text x="45" y="126" class="big">{percent:.1f}%</text>
<text x="45" y="160" class="statement">A majority of Culsans neighbourhoods contain predicted</text>
<text x="45" y="190" class="statement">defence or anti-defence systems</text>
<text x="48" y="217" class="note">{associated} of {total} analysable representative RefSeq neighbourhoods</text>
<line x1="45" y1="240" x2="955" y2="240" stroke="#D5DCE5" stroke-width="2"/>
<text x="70" y="291" class="metric">{len(genera)}</text><text x="70" y="318" class="label">bacterial genera represented</text>
<text x="370" y="291" class="metric">{multiple}</text><text x="370" y="318" class="label">positive neighbourhoods with ≥2 nonredundant loci</text>
<text x="740" y="291" class="metric">{multiple_percent:.1f}%</text><text x="740" y="318" class="label">of positive neighbourhoods</text>
</svg>'''


def make_environment_percentage_svg(associated: int, total: int, colour: str) -> str:
    percent = 100 * associated / total if total else 0.0
    remainder = total - associated
    width = 900
    bar_x, bar_y, bar_width, bar_height = 70, 115, 760, 62
    yes_width = bar_width * percent / 100
    no_colour = tint(colour, 0.84)
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="280" viewBox="0 0 {width} 280">
<style>.title{{font:700 25px Arial,sans-serif;fill:#111}}.note{{font:14px Arial,sans-serif;fill:#444}}.inside{{font:700 18px Arial,sans-serif}}.legend{{font:15px Arial,sans-serif;fill:#111}}</style>
<rect width="100%" height="100%" fill="#FFFFFF"/>
<text x="70" y="42" class="title">Culsans in predicted defence environments</text>
<text x="70" y="70" class="note">One representative RefSeq neighbourhood per analysable protein accession (n={total})</text>
<rect x="{bar_x}" y="{bar_y}" width="{bar_width}" height="{bar_height}" fill="{no_colour}" stroke="#111" stroke-width="2"/>
<rect x="{bar_x}" y="{bar_y}" width="{yes_width:.2f}" height="{bar_height}" fill="{colour}" stroke="#111" stroke-width="2"/>
<text x="{bar_x + yes_width / 2:.2f}" y="{bar_y + 39}" text-anchor="middle" class="inside" fill="#FFFFFF">{associated} ({percent:.2f}%)</text>
<text x="{bar_x + yes_width + (bar_width - yes_width) / 2:.2f}" y="{bar_y + 39}" text-anchor="middle" class="inside" fill="#111111">{remainder} ({100 - percent:.2f}%)</text>
<rect x="70" y="215" width="18" height="18" fill="{colour}" stroke="#111"/><text x="98" y="229" class="legend">At least one defence/anti-defence system</text>
<rect x="480" y="215" width="18" height="18" fill="{no_colour}" stroke="#111"/><text x="508" y="229" class="legend">No system detected</text>
</svg>''']
    return "\n".join(parts)


def joined_unique(values: Iterable[object]) -> str:
    return "; ".join(dict.fromkeys(str(value).strip() for value in values if str(value).strip()))


def make_supplementary_tables(
    combined: list[dict[str, str]],
    contexts: list[dict[str, object]],
    summary_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Combine scanner calls with gene maps and genomic coordinates."""
    context_by_protein = {str(row["protein_accession"]): row for row in contexts}
    summary_by_protein = {str(row["protein_accession"]): row for row in summary_rows}
    calls_by_protein: dict[str, list[dict[str, str]]] = defaultdict(list)
    detail_rows: list[dict[str, object]] = []
    gene_cache: dict[str, tuple[dict[str, dict[str, str]], dict[str, str] | None]] = {}

    for call in combined:
        protein = call["protein_accession"]
        calls_by_protein[protein].append(call)
        context = context_by_protein[protein]
        if protein not in gene_cache:
            gene_rows = read_delimited(Path(str(context["context_dir"])) / "genes.tsv", "\t")
            gene_map = {row["feature_id"]: row for row in gene_rows}
            target = next((row for row in gene_rows if row.get("is_target") == "yes"), None)
            gene_cache[protein] = (gene_map, target)
        gene_map, target = gene_cache[protein]
        raw = json.loads(call["raw_result"])
        if call["tool"] == "PADLOC":
            members = raw if isinstance(raw, list) else [raw]
            feature_ids = list(dict.fromkeys(
                member.get("target.name", "") for member in members if member.get("target.name")
            ))
            profiles = joined_unique(member.get("hmm.name", "") for member in members)
            accessions = joined_unique(member.get("hmm.accession", "") for member in members)
            reported_gene_count = len(members)
        else:
            members = [raw]
            feature_ids = list(dict.fromkeys(re.findall(r"FG_\d+", str(raw.get("protein_in_syst", "")))))
            profiles = str(raw.get("name_of_profiles_in_sys", ""))
            accessions = ""
            reported_gene_count = int(raw.get("genes_count", 0) or 0)
        mapped = [gene_map[feature_id] for feature_id in feature_ids if feature_id in gene_map]
        missing_ids = [feature_id for feature_id in feature_ids if feature_id not in gene_map]
        local_start = min((int(gene["local_start"]) for gene in mapped), default=None)
        local_stop = max((int(gene["local_stop"]) for gene in mapped), default=None)
        target_start = int(target["local_start"]) if target else None
        target_stop = int(target["local_stop"]) if target else None
        target_strand = target["strand"] if target else ""
        gap: int | None = None
        relative = "unavailable"
        if local_start is not None and local_stop is not None and target_start is not None and target_stop is not None:
            if local_stop < target_start:
                gap = target_start - local_stop - 1
                relative = "upstream" if target_strand == "+" else "downstream"
            elif local_start > target_stop:
                gap = local_start - target_stop - 1
                relative = "downstream" if target_strand == "+" else "upstream"
            else:
                gap = 0
                relative = "overlaps Culsans CDS"
        gene_local_coordinates = joined_unique(
            f"{gene['feature_id']}:{gene['local_start']}-{gene['local_stop']}({gene['strand']})" for gene in mapped
        )
        gene_genomic_coordinates = joined_unique(
            f"{gene['feature_id']}:{gene['genomic_start']}-{gene['genomic_stop']}({gene['strand']})" for gene in mapped
        )
        detail_rows.append({
            "Culsans_protein_accession": protein,
            "RefSeq_nucleotide_accession": call["nuccore_accession"],
            "organism": context["organism"],
            "strain": context["strain"],
            "assembly": context["assembly"],
            "program": call["tool"],
            "call_class": "anti-defence" if call["tool"] == "AntiDefenseFinder" else "defence",
            "system": call["system"],
            "system_id": call["system_id"],
            "model_profiles": profiles,
            "model_accessions": accessions,
            "reported_gene_count": reported_gene_count,
            "mapped_gene_count": len(mapped),
            "local_feature_ids": joined_unique(feature_ids),
            "RefSeq_system_protein_accessions": joined_unique(gene["original_protein_id"] for gene in mapped),
            "locus_tags": joined_unique(gene["locus_tag"] for gene in mapped),
            "gene_products": joined_unique(gene["product"] for gene in mapped),
            "gene_strands": joined_unique(gene["strand"] for gene in mapped),
            "gene_coordinates_in_neighborhood": gene_local_coordinates,
            "gene_coordinates_on_replicon": gene_genomic_coordinates,
            "system_local_start": local_start if local_start is not None else "",
            "system_local_stop": local_stop if local_stop is not None else "",
            "position_relative_to_Culsans": relative,
            "nearest_distance_to_Culsans_bp": gap if gap is not None else "",
            "Culsans_local_start": target_start if target_start is not None else "",
            "Culsans_local_stop": target_stop if target_stop is not None else "",
            "Culsans_genomic_start": target["genomic_start"] if target else "",
            "Culsans_genomic_stop": target["genomic_stop"] if target else "",
            "Culsans_strand": target_strand,
            "retrieved_intervals": context["request_intervals"],
            "replicon_topology": context["replicon_topology"],
            "edge_status": context["edge_status"],
            "coordinate_mapping_status": "complete" if mapped and not missing_ids else (
                "partial: " + joined_unique(missing_ids) if mapped else "unmapped: " + joined_unique(missing_ids)
            ),
        })

    context_rows: list[dict[str, object]] = []
    for protein, context in context_by_protein.items():
        summary = summary_by_protein[protein]
        calls = calls_by_protein.get(protein, [])
        by_tool = {tool: [row for row in calls if row["tool"] == tool] for tool in (
            "DefenseFinder", "AntiDefenseFinder", "PADLOC",
        )}
        context_rows.append({
            "Culsans_protein_accession": protein,
            "RefSeq_nucleotide_accession": context["nuccore_accession"],
            "organism": context["organism"],
            "strain": context["strain"],
            "assembly": context["assembly"],
            "Culsans_CDS_start": context["cds_start"],
            "Culsans_CDS_stop": context["cds_stop"],
            "retrieved_intervals": context["request_intervals"],
            "replicon_topology": context["replicon_topology"],
            "edge_status": context["edge_status"],
            "defence_or_antidefence_detected": summary["defence_environment"],
            "DefenseFinder_call_count": len(by_tool["DefenseFinder"]),
            "DefenseFinder_systems": joined_unique(row["system"] for row in by_tool["DefenseFinder"]),
            "AntiDefenseFinder_call_count": len(by_tool["AntiDefenseFinder"]),
            "AntiDefenseFinder_systems": joined_unique(row["system"] for row in by_tool["AntiDefenseFinder"]),
            "PADLOC_call_count": len(by_tool["PADLOC"]),
            "PADLOC_systems": joined_unique(row["system"] for row in by_tool["PADLOC"]),
            "combined_tool_specific_call_count": len(calls),
            "combined_tool_specific_systems": joined_unique(f"{row['tool']}: {row['system']}" for row in calls),
            "target_Culsans_CDS_found": summary["target_cds_found"],
        })
    return detail_rows, context_rows


def normalized_system_name(name: str) -> str:
    value = name.lower().replace("_", "-")
    value = re.sub(r"[^a-z0-9]+", "-", value).strip("-")
    value = re.sub(
        r"-type-(iis|iv|iii|ii|i)(?=-|$)",
        lambda match: "-" + {"iis": "ii", "iv": "iv", "iii": "iii", "ii": "ii", "i": "i"}[match.group(1)],
        value,
    )
    aliases = {
        "septu": "septu-i",
        "cas-class1-subtype-i-c": "cas-i-c",
        "gao-hhe": "hhe",
        "gao-ppl": "ppl",
        "pago-s1b": "argonaute-i",
    }
    return aliases.get(value, value)


def canonical_system_name(name: str) -> str:
    """Return one stable biological label for equivalent program nomenclature."""
    key = normalized_system_name(name)
    explicit = {
        "cbass-i": "CBASS type I", "cbass-ii": "CBASS type II",
        "rm-i": "RM type I", "rm-ii": "RM type II", "rm-iii": "RM type III",
        "rm-iv": "RM type IV", "rm-iv-1": "RM type IV-1",
        "gabija": "Gabija", "septu-i": "Septu type I",
        "cas-i-c": "CRISPR-Cas type I-C", "wadjet-i": "Wadjet type I",
        "hhe": "Hhe", "ppl": "Ppl", "argonaute-i": "Argonaute type I",
    }
    if key in explicit:
        return explicit[key]
    return name.replace("_", " ")


def preferred_system_name(rows: list[dict[str, object]]) -> str:
    names_by_key: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in rows:
        names_by_key[normalized_system_name(str(row["system"]))].append((str(row["program"]), str(row["system"])))
    displayed = []
    for key, values in sorted(names_by_key.items()):
        ordered = sorted(values, key=lambda pair: ({"DefenseFinder": 0, "AntiDefenseFinder": 1, "PADLOC": 2}.get(pair[0], 9), pair[1]))
        name = canonical_system_name(ordered[0][1])
        displayed.append(name)
    return " / ".join(displayed)


def make_consensus_loci(
    system_calls: list[dict[str, object]],
    total_neighborhoods: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Merge overlapping tool-specific calls into nonredundant genomic loci."""
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in system_calls:
        if normalized_system_name(str(row["system"])) == "dms-other":
            continue
        grouped[(str(row["Culsans_protein_accession"]), str(row["call_class"]))].append(row)
    consensus: list[dict[str, object]] = []
    for (protein, call_class), rows in grouped.items():
        parents = list(range(len(rows)))

        def find(index: int) -> int:
            while parents[index] != index:
                parents[index] = parents[parents[index]]
                index = parents[index]
            return index

        def union(left: int, right: int) -> None:
            root_left, root_right = find(left), find(right)
            if root_left != root_right:
                parents[root_right] = root_left

        feature_sets = [set(str(row["local_feature_ids"]).split("; ")) - {""} for row in rows]
        for left in range(len(rows)):
            for right in range(left + 1, len(rows)):
                same_family = normalized_system_name(str(rows[left]["system"])) == normalized_system_name(str(rows[right]["system"]))
                if same_family and feature_sets[left] & feature_sets[right]:
                    union(left, right)
                    continue
                # Overlap alone is not enough: unrelated models can share a
                # component gene and must remain distinct biological calls.
                left_start, left_stop = rows[left]["system_local_start"], rows[left]["system_local_stop"]
                right_start, right_stop = rows[right]["system_local_start"], rows[right]["system_local_stop"]
                if same_family and all(value != "" for value in (left_start, left_stop, right_start, right_stop)):
                    if max(int(left_start), int(right_start)) <= min(int(left_stop), int(right_stop)):
                        union(left, right)
        clusters: dict[int, list[dict[str, object]]] = defaultdict(list)
        for index, row in enumerate(rows):
            clusters[find(index)].append(row)
        ordered_clusters = sorted(clusters.values(), key=lambda cluster: (
            min(int(row["system_local_start"]) for row in cluster if row["system_local_start"] != ""),
            preferred_system_name(cluster),
        ))
        for locus_number, cluster in enumerate(ordered_clusters, 1):
            programs = sorted({str(row["program"]) for row in cluster})
            local_starts = [int(row["system_local_start"]) for row in cluster if row["system_local_start"] != ""]
            local_stops = [int(row["system_local_stop"]) for row in cluster if row["system_local_stop"] != ""]
            distances = [int(row["nearest_distance_to_Culsans_bp"]) for row in cluster if row["nearest_distance_to_Culsans_bp"] != ""]
            first = cluster[0]
            consensus.append({
                "Culsans_protein_accession": protein,
                "RefSeq_nucleotide_accession": first["RefSeq_nucleotide_accession"],
                "organism": first["organism"],
                "strain": first["strain"],
                "assembly": first["assembly"],
                "consensus_locus_id": f"{protein}_L{locus_number:03d}",
                "call_class": call_class,
                "consensus_system": preferred_system_name(cluster),
                "supporting_programs": joined_unique(programs),
                "program_agreement": "DefenseFinder and PADLOC" if {"DefenseFinder", "PADLOC"}.issubset(programs) else joined_unique(programs),
                "DefenseFinder_systems": joined_unique(row["system"] for row in cluster if row["program"] == "DefenseFinder"),
                "AntiDefenseFinder_systems": joined_unique(row["system"] for row in cluster if row["program"] == "AntiDefenseFinder"),
                "PADLOC_systems": joined_unique(row["system"] for row in cluster if row["program"] == "PADLOC"),
                "tool_specific_call_count": len(cluster),
                "local_feature_ids": joined_unique(value for row in cluster for value in str(row["local_feature_ids"]).split("; ")),
                "RefSeq_system_protein_accessions": joined_unique(value for row in cluster for value in str(row["RefSeq_system_protein_accessions"]).split("; ")),
                "locus_tags": joined_unique(value for row in cluster for value in str(row["locus_tags"]).split("; ")),
                "gene_products": joined_unique(value for row in cluster for value in str(row["gene_products"]).split("; ")),
                "gene_coordinates_in_neighborhood": joined_unique(value for row in cluster for value in str(row["gene_coordinates_in_neighborhood"]).split("; ")),
                "gene_coordinates_on_replicon": joined_unique(value for row in cluster for value in str(row["gene_coordinates_on_replicon"]).split("; ")),
                "consensus_local_start": min(local_starts) if local_starts else "",
                "consensus_local_stop": max(local_stops) if local_stops else "",
                "position_relative_to_Culsans": joined_unique(row["position_relative_to_Culsans"] for row in cluster),
                "nearest_distance_to_Culsans_bp": min(distances) if distances else "",
                "Culsans_genomic_start": first["Culsans_genomic_start"],
                "Culsans_genomic_stop": first["Culsans_genomic_stop"],
                "Culsans_strand": first["Culsans_strand"],
                "retrieved_intervals": first["retrieved_intervals"],
                "edge_status": first["edge_status"],
            })
    frequency_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in consensus:
        frequency_sets[(str(row["call_class"]), str(row["consensus_system"]))].add(str(row["Culsans_protein_accession"]))
    frequency_rows = [{
        "call_class": call_class,
        "consensus_system": system,
        "neighborhood_count": len(proteins),
        "percent_of_analysed_neighborhoods": f"{100 * len(proteins) / total_neighborhoods:.2f}" if total_neighborhoods else "0.00",
    } for (call_class, system), proteins in frequency_sets.items()]
    frequency_rows.sort(key=lambda row: (-int(row["neighborhood_count"]), str(row["call_class"]), str(row["consensus_system"])))
    return consensus, frequency_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hit_table", type=Path, help="BLAST hit table containing RefSeq protein accessions")
    parser.add_argument("outdir", type=Path, help="Output directory")
    parser.add_argument("--email", required=True, help="Contact email required by NCBI E-utilities")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"), help="NCBI API key or NCBI_API_KEY")
    parser.add_argument("--protein-column", help="Header of the protein-accession column; detected automatically by default")
    parser.add_argument("--flank-bp", type=int, default=30000, help="Bases extracted on each side of the target CDS (default: 30000)")
    parser.add_argument("--defensefinder-env", default="defensefinder", help="Conda environment containing DefenseFinder")
    parser.add_argument("--padloc-env", default="padloc", help="Conda environment containing PADLOC")
    parser.add_argument("--conda", default=shutil.which("conda") or shutil.which("mamba") or "conda", help="Conda or Mamba executable")
    parser.add_argument("--threads", type=int, default=1, help="CPUs passed to each scanner (default: 1)")
    parser.add_argument("--jobs", type=int, default=1, help="Neighborhoods scanned concurrently (default: 1)")
    parser.add_argument("--heatmap-color", default="#264979", help="Six-digit hex colour used for plots (default: #264979)")
    parser.add_argument("--heatmap-top-systems", type=int, default=25, help="Maximum system columns in the heatmap")
    parser.add_argument("--heatmap-max-neighborhoods", type=int, default=250, help="Maximum neighbourhood rows in the heatmap")
    parser.add_argument("--max-proteins", type=int, help="Limit unique proteins for a test run")
    parser.add_argument("--delay", type=float, default=None, help="Seconds between NCBI requests")
    parser.add_argument("--offline", action="store_true", help="Use only cached NCBI responses")
    parser.add_argument("--dry-run", action="store_true", help="Prepare inputs and print scanner commands without running them")
    parser.add_argument("--skip-tool-checks", action="store_true", help="Skip scanner environment preflight checks")
    parser.add_argument("--skip-defensefinder", action="store_true", help="Do not run or parse DefenseFinder")
    parser.add_argument("--skip-padloc", action="store_true", help="Do not run or parse PADLOC")
    args = parser.parse_args()

    if not args.hit_table.is_file():
        parser.error(f"Hit table not found: {args.hit_table}")
    if args.flank_bp < 0 or args.threads < 1 or args.jobs < 1:
        parser.error("--flank-bp must be non-negative; --threads and --jobs must be at least 1")
    try:
        hex_to_rgb(args.heatmap_color)
    except ValueError as exc:
        parser.error(str(exc))
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not args.dry_run and not args.skip_tool_checks:
        print("Checking scanner environments...", flush=True)
        try:
            if not args.skip_defensefinder:
                check_tool(
                    conda_prefix(args.conda, args.defensefinder_env) + ["defense-finder", "run", "--help"],
                    "DefenseFinder",
                )
            if not args.skip_padloc:
                check_tool(
                    conda_prefix(args.conda, args.padloc_env) + ["padloc", "--check-deps"],
                    "PADLOC", reject_missing_warning=True,
                )
        except RuntimeError as exc:
            parser.error(str(exc))
    cache = args.outdir / "cache"
    delay = args.delay if args.delay is not None else (0.11 if args.api_key else 0.34)
    client = NcbiClient(args.email, args.api_key, delay)
    proteins, input_rows, protein_column = read_proteins(args.hit_table, args.protein_column)
    if args.max_proteins:
        proteins = proteins[:args.max_proteins]
    print(f"Input rows: {input_rows}; unique proteins selected: {len(proteins)}; column: {protein_column}")

    representatives: list[dict[str, str]] = []
    failures: list[dict[str, str]] = []
    for index, protein in enumerate(proteins, 1):
        print(f"[IPG {index}/{len(proteins)}] {protein}", flush=True)
        try:
            xml = cache_fetch(client, cache / "ipg" / f"{safe_name(protein)}.xml", "efetch.fcgi", {
                "db": "protein", "id": protein, "rettype": "ipg", "retmode": "xml",
            }, args.offline)
            representatives.append(choose_representative(parse_ipg_xml(xml, protein)))
        except Exception as exc:
            print(f"  FAILED: {exc}", file=sys.stderr)
            failures.append({"stage": "ipg", "protein_accession": protein, "error": str(exc)})

    representative_fields = [
        "protein_accession", "ipg_protein_accession", "protein_name", "nuccore_accession",
        "cds_start", "cds_stop", "strand", "organism", "strain", "assembly",
        "ipg_refseq_candidates", "representative_rank", "representative_rule",
    ]
    write_tsv(args.outdir / "representatives.tsv", representatives, representative_fields)

    contexts: list[dict[str, object]] = []
    for index, representative in enumerate(representatives, 1):
        protein = representative["protein_accession"]
        accession = representative["nuccore_accession"]
        print(f"[Context {index}/{len(representatives)}] {protein} -> {accession}", flush=True)
        context_dir = args.outdir / "neighborhoods" / safe_name(protein)
        try:
            returned_accession, sequence, features, metadata = fetch_context(
                client, representative, args.flank_bp, cache, args.offline,
            )
            write_context_files(context_dir, representative, returned_accession, sequence, features)
            context_row: dict[str, object] = {**representative, **metadata, "context_dir": str(context_dir)}
            contexts.append(context_row)
            atomic_write(context_dir / "context.json", json.dumps(context_row, indent=2) + "\n")
        except Exception as exc:
            print(f"  FAILED: {exc}", file=sys.stderr)
            failures.append({"stage": "context", "protein_accession": protein, "error": str(exc)})

    context_fields = representative_fields + [
        "returned_accession", "definition", "organism_from_record", "request_start",
        "requested_stop", "request_intervals", "returned_length_bp", "returned_stop",
        "coordinate_mode", "target_cds_found", "cds_count", "replicon_length",
        "replicon_topology", "edge_status", "context_dir",
    ]
    write_tsv(args.outdir / "neighborhood_manifest.tsv", contexts, context_fields)

    def scan_context(context: dict[str, object]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
        found: list[dict[str, str]] = []
        scan_failures: list[dict[str, str]] = []
        protein = str(context["protein_accession"])
        context_dir = Path(str(context["context_dir"]))
        faa = context_dir / "neighborhood.faa"
        gff = context_dir / "neighborhood.gff"
        if not args.skip_defensefinder:
            df_dir = context_dir / "defensefinder"
            marker = df_dir / ".flankguard_complete"
            try:
                if not marker.exists() and not args.dry_run:
                    df_dir.mkdir(parents=True, exist_ok=True)
                if not marker.exists():
                    command = conda_prefix(args.conda, args.defensefinder_env) + [
                        "defense-finder", "run", str(faa), "--out-dir", str(df_dir),
                        "--workers", str(args.threads), "--antidefensefinder",
                        "--index-dir", str(df_dir / "index"),
                    ]
                    run_logged(command, df_dir / "flankguard.log", args.dry_run)
                    if not args.dry_run:
                        atomic_write(marker, "complete\n")
                if not args.dry_run:
                    for result in parse_defensefinder(df_dir):
                        found.append({
                            "protein_accession": protein,
                            "nuccore_accession": str(context["nuccore_accession"]),
                            **result,
                        })
            except Exception as exc:
                print(f"  DefenseFinder FAILED: {exc}", file=sys.stderr)
                scan_failures.append({"stage": "defensefinder", "protein_accession": protein, "error": str(exc)})
        if not args.skip_padloc:
            padloc_dir = context_dir / "padloc"
            marker = padloc_dir / ".flankguard_complete"
            try:
                if not marker.exists() and not args.dry_run:
                    padloc_dir.mkdir(parents=True, exist_ok=True)
                if not marker.exists():
                    command = conda_prefix(args.conda, args.padloc_env) + [
                        "padloc", "--faa", str(faa), "--gff", str(gff),
                        "--outdir", str(padloc_dir), "--cpu", str(args.threads),
                    ]
                    run_logged(command, padloc_dir / "flankguard.log", args.dry_run)
                    if not args.dry_run:
                        atomic_write(marker, "complete\n")
                if not args.dry_run:
                    for result in parse_padloc(padloc_dir):
                        found.append({
                            "protein_accession": protein,
                            "nuccore_accession": str(context["nuccore_accession"]),
                            **result,
                        })
            except Exception as exc:
                print(f"  PADLOC FAILED: {exc}", file=sys.stderr)
                scan_failures.append({"stage": "padloc", "protein_accession": protein, "error": str(exc)})
        return found, scan_failures

    combined: list[dict[str, str]] = []
    if args.jobs == 1:
        for index, context in enumerate(contexts, 1):
            print(f"[Scan {index}/{len(contexts)}] {context['protein_accession']}", flush=True)
            found, scan_failures = scan_context(context)
            combined.extend(found)
            failures.extend(scan_failures)
    else:
        print(f"Scanning {len(contexts)} neighborhoods with {args.jobs} concurrent jobs", flush=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            future_to_context = {executor.submit(scan_context, context): context for context in contexts}
            completed_count = 0
            for future in as_completed(future_to_context):
                completed_count += 1
                context = future_to_context[future]
                protein = str(context["protein_accession"])
                try:
                    found, scan_failures = future.result()
                    combined.extend(found)
                    failures.extend(scan_failures)
                except Exception as exc:
                    failures.append({"stage": "scan", "protein_accession": protein, "error": str(exc)})
                print(f"[Scan {completed_count}/{len(contexts)} complete] {protein}", flush=True)

    combined_fields = [
        "protein_accession", "nuccore_accession", "tool", "system", "system_id",
        "result_file", "raw_result",
    ]
    write_tsv(args.outdir / "combined_systems.tsv", combined, combined_fields)
    analysis_combined = [
        row for row in combined
        if normalized_system_name(str(row["system"])) != "dms-other"
    ]
    presence: dict[str, set[str]] = defaultdict(set)
    tools_by_protein: dict[str, set[str]] = defaultdict(set)
    systems_by_protein: dict[str, set[str]] = defaultdict(set)
    for row in analysis_combined:
        label = f"{row['tool']}: {row['system']}"
        presence[row["protein_accession"]].add(label)
        tools_by_protein[row["protein_accession"]].add(row["tool"])
        systems_by_protein[row["protein_accession"]].add(label)

    summary_rows: list[dict[str, object]] = []
    for context in contexts:
        protein = str(context["protein_accession"])
        labels = sorted(systems_by_protein.get(protein, set()))
        tools = sorted(tools_by_protein.get(protein, set()))
        summary_rows.append({
            "protein_accession": protein,
            "nuccore_accession": context["nuccore_accession"],
            "organism": context["organism"],
            "assembly": context["assembly"],
            "defence_environment": "yes" if labels else "no",
            "total_systems": len(labels),
            "tools_with_calls": "; ".join(tools),
            "systems": "; ".join(labels),
            "target_cds_found": context["target_cds_found"],
        })
    summary_fields = [
        "protein_accession", "nuccore_accession", "organism", "assembly",
        "defence_environment", "total_systems", "tools_with_calls", "systems",
        "target_cds_found",
    ]
    write_tsv(args.outdir / "neighborhood_summary.tsv", summary_rows, summary_fields)

    supplementary_calls, supplementary_contexts = make_supplementary_tables(combined, contexts, summary_rows)
    consensus_loci, consensus_frequencies = make_consensus_loci(supplementary_calls, len(summary_rows))
    consensus_by_protein: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in consensus_loci:
        consensus_by_protein[str(row["Culsans_protein_accession"])].append(row)
    for row in supplementary_contexts:
        loci = consensus_by_protein.get(str(row["Culsans_protein_accession"]), [])
        row["nonredundant_consensus_locus_count"] = len(loci)
        row["nonredundant_consensus_systems"] = joined_unique(locus["consensus_system"] for locus in loci)
    supplementary_call_fields = [
        "Culsans_protein_accession", "RefSeq_nucleotide_accession", "organism", "strain", "assembly",
        "program", "call_class", "system", "system_id", "model_profiles", "model_accessions",
        "reported_gene_count", "mapped_gene_count", "local_feature_ids",
        "RefSeq_system_protein_accessions", "locus_tags", "gene_products", "gene_strands",
        "gene_coordinates_in_neighborhood", "gene_coordinates_on_replicon",
        "system_local_start", "system_local_stop", "position_relative_to_Culsans",
        "nearest_distance_to_Culsans_bp", "Culsans_local_start", "Culsans_local_stop",
        "Culsans_genomic_start", "Culsans_genomic_stop", "Culsans_strand",
        "retrieved_intervals", "replicon_topology", "edge_status", "coordinate_mapping_status",
    ]
    supplementary_context_fields = [
        "Culsans_protein_accession", "RefSeq_nucleotide_accession", "organism", "strain", "assembly",
        "Culsans_CDS_start", "Culsans_CDS_stop", "retrieved_intervals", "replicon_topology", "edge_status",
        "defence_or_antidefence_detected", "DefenseFinder_call_count", "DefenseFinder_systems",
        "AntiDefenseFinder_call_count", "AntiDefenseFinder_systems", "PADLOC_call_count", "PADLOC_systems",
        "combined_tool_specific_call_count", "combined_tool_specific_systems", "target_Culsans_CDS_found",
        "nonredundant_consensus_locus_count", "nonredundant_consensus_systems",
    ]
    consensus_locus_fields = [
        "Culsans_protein_accession", "RefSeq_nucleotide_accession", "organism", "strain", "assembly",
        "consensus_locus_id", "call_class", "consensus_system", "supporting_programs", "program_agreement",
        "DefenseFinder_systems", "AntiDefenseFinder_systems", "PADLOC_systems", "tool_specific_call_count",
        "local_feature_ids", "RefSeq_system_protein_accessions", "locus_tags", "gene_products",
        "gene_coordinates_in_neighborhood", "gene_coordinates_on_replicon", "consensus_local_start",
        "consensus_local_stop", "position_relative_to_Culsans", "nearest_distance_to_Culsans_bp",
        "Culsans_genomic_start", "Culsans_genomic_stop", "Culsans_strand", "retrieved_intervals", "edge_status",
    ]
    write_tsv(args.outdir / "supplementary_system_calls.tsv", supplementary_calls, supplementary_call_fields)
    write_tsv(args.outdir / "supplementary_neighborhoods.tsv", supplementary_contexts, supplementary_context_fields)
    write_tsv(args.outdir / "supplementary_consensus_loci.tsv", consensus_loci, consensus_locus_fields)
    write_tsv(args.outdir / "consensus_system_frequencies.tsv", consensus_frequencies, [
        "call_class", "consensus_system", "neighborhood_count", "percent_of_analysed_neighborhoods",
    ])

    total = len(summary_rows)
    associated = sum(row["defence_environment"] == "yes" for row in summary_rows)
    environment_rows = [
        {
            "category": "At least one detected defence/anti-defence system",
            "count": associated,
            "percent": f"{100 * associated / total:.2f}" if total else "0.00",
        },
        {
            "category": "No system detected",
            "count": total - associated,
            "percent": f"{100 * (total - associated) / total:.2f}" if total else "0.00",
        },
    ]
    write_tsv(args.outdir / "query_environment_counts.tsv", environment_rows, ["category", "count", "percent"])

    frequency_sets: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in combined:
        frequency_sets[(row["tool"], row["system"])].add(row["protein_accession"])
    frequency_rows = [
        {
            "tool": tool,
            "system": system,
            "neighborhood_count": len(protein_set),
            "percent_of_neighborhoods": f"{100 * len(protein_set) / total:.2f}" if total else "0.00",
        }
        for (tool, system), protein_set in frequency_sets.items()
    ]
    frequency_rows.sort(key=lambda row: (-int(row["neighborhood_count"]), str(row["tool"]), str(row["system"])))
    write_tsv(args.outdir / "system_frequencies.tsv", frequency_rows, [
        "tool", "system", "neighborhood_count", "percent_of_neighborhoods",
    ])
    tool_frequencies = Counter({f"{row['tool']}: {row['system']}": int(row["neighborhood_count"]) for row in frequency_rows})
    consensus_presence: dict[str, set[str]] = defaultdict(set)
    for locus in consensus_loci:
        consensus_presence[str(locus["Culsans_protein_accession"])].add(str(locus["consensus_system"]))
    consensus_frequency_counts = Counter({
        str(row["consensus_system"]): int(row["neighborhood_count"])
        for row in consensus_frequencies
    })
    top_labels = [label for label, _ in consensus_frequency_counts.most_common(args.heatmap_top_systems)]
    matrix_rows: list[dict[str, object]] = []
    summary_by_protein = {str(row["protein_accession"]): row for row in summary_rows}
    for protein in sorted(summary_by_protein):
        summary = summary_by_protein[protein]
        matrix_rows.append({
            "protein_accession": protein,
            "nuccore_accession": summary["nuccore_accession"],
            "organism": summary["organism"],
            **{label: int(label in consensus_presence.get(protein, set())) for label in top_labels},
        })
    write_tsv(args.outdir / "system_presence_matrix.tsv", matrix_rows, [
        "protein_accession", "nuccore_accession", "organism", *top_labels,
    ])
    atomic_write(args.outdir / "system_heatmap.svg", make_heatmap(
        summary_rows, consensus_presence, consensus_frequency_counts, args.heatmap_color,
        args.heatmap_top_systems, args.heatmap_max_neighborhoods,
    ))
    atomic_write(args.outdir / "tool_specific_system_heatmap.svg", make_heatmap(
        summary_rows, presence, tool_frequencies, args.heatmap_color,
        args.heatmap_top_systems, args.heatmap_max_neighborhoods,
    ))
    atomic_write(args.outdir / "system_heatmap_A4.svg", make_a4_taxon_heatmap(
        summary_rows, consensus_presence, consensus_frequency_counts, args.heatmap_color,
    ))
    atomic_write(args.outdir / "tool_specific_system_frequency_strip.svg", make_frequency_heatmap_svg(
        frequency_rows, args.heatmap_color, args.heatmap_top_systems,
    ))
    atomic_write(args.outdir / "tool_specific_system_frequencies.svg", make_frequency_svg(
        frequency_rows, args.heatmap_color, args.heatmap_top_systems,
    ))
    atomic_write(args.outdir / "system_frequencies.svg", make_frequency_svg(
        consensus_frequencies, args.heatmap_color, args.heatmap_top_systems, tool_specific=False,
    ))
    atomic_write(args.outdir / "defence_environment_percentage.svg", make_environment_percentage_svg(
        associated, total, args.heatmap_color,
    ))
    atomic_write(args.outdir / "association_breadth.svg", make_association_breadth_svg(
        summary_rows, consensus_loci, args.heatmap_color,
    ))
    write_tsv(args.outdir / "failures.tsv", failures, ["stage", "protein_accession", "error"])
    run_summary = {
        "input": str(args.hit_table.resolve()),
        "input_rows": input_rows,
        "protein_column": protein_column,
        "unique_proteins_selected": len(proteins),
        "representatives_selected": len(representatives),
        "neighborhoods_extracted": len(contexts),
        "flank_bp_each_side": args.flank_bp,
        "neighborhoods_with_detected_systems": associated,
        "neighborhoods_without_detected_systems": total - associated,
        "defensefinder_environment": args.defensefinder_env,
        "padloc_environment": args.padloc_env,
        "heatmap_color": args.heatmap_color,
        "dry_run": args.dry_run,
        "failures": len(failures),
    }
    atomic_write(args.outdir / "summary.json", json.dumps(run_summary, indent=2) + "\n")
    print(f"Done. Combined table: {args.outdir / 'neighborhood_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
