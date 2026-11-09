#!/usr/bin/env python3

"""
Clean the genome set for the PAPI29/Culsans Pseudomonas analysis.

INPUTS
------
genomes.txt
    Current genome list. This contains:
      - the original dereplicated Pseudomonas background genomes
      - Culsans representative genomes that were subsequently added

papi29distri
    IPG table containing Culsans protein hits, with columns including:
      blast_protein_accession
      ipg_source
      ipg_protein_accession
      protein_name
      nuccore_accession
      cds_start
      cds_stop
      strand
      organism
      strain
      assembly

ncbi/
    Directory containing downloaded genome files.

OUTPUTS
-------
genomes_pseudomonas.txt
    Clean genome list for downstream GToTree analysis.

genomes_pseudomonas_report.tsv
    Audit trail showing removed, retained, added and missing
    PAPI29 assemblies.
"""

from pathlib import Path
import csv
import re
from collections import defaultdict


# ============================================================
# SETTINGS
# ============================================================

GENOMES_FILE = Path("genomes.txt")
IPG_FILE = Path("papi29distri")
NCBI_DIR = Path("ncbi")

OUTPUT = Path("genomes_pseudomonas.txt")
REPORT = Path("genomes_pseudomonas_report.tsv")

# Match things such as:
# GCF_000284555.1
# GCA_012345678.2
ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


# ============================================================
# FUNCTIONS
# ============================================================

def get_accession(text):
    """Extract GCF/GCA assembly accession from a path/string."""
    match = ACC_RE.search(str(text))

    if match:
        return match.group(1)

    return None


def is_pseudomonas(organism):
    """
    Return True only when organism explicitly starts with Pseudomonas.
    Case-insensitive.
    """
    organism = organism.strip().lower()

    return organism.startswith("pseudomonas ")


# ============================================================
# CHECK INPUTS
# ============================================================

for required in [GENOMES_FILE, IPG_FILE, NCBI_DIR]:

    if not required.exists():
        raise SystemExit(
            f"ERROR: required input does not exist: {required}"
        )


# ============================================================
# 1. READ THE IPG FILE
# ============================================================

# We retain all organism assignments for an assembly rather than
# immediately assuming one row = one assembly classification.

ipg_organisms = defaultdict(set)

with IPG_FILE.open(newline="") as handle:

    reader = csv.DictReader(handle, delimiter="\t")

    required_columns = {"assembly", "organism"}

    if not required_columns.issubset(reader.fieldnames or []):
        raise SystemExit(
            "ERROR: papi29distri does not contain the required "
            "'assembly' and 'organism' columns.\n"
            f"Columns found: {reader.fieldnames}"
        )

    for row in reader:

        assembly = row["assembly"].strip()
        organism = row["organism"].strip()

        # Skip IPG rows without a usable assembly accession.
        if not assembly:
            continue

        # Extract GCF/GCA in case anything extra appears in field.
        assembly_acc = get_accession(assembly)

        if assembly_acc is None:
            continue

        if organism:
            ipg_organisms[assembly_acc].add(organism)


# ============================================================
# 2. CLASSIFY PAPI29 ASSEMBLIES FROM THE IPG
# ============================================================

papi29_pseudomonas = {}
papi29_nonpseudomonas = {}
ambiguous_ipg = {}


for assembly, organisms in ipg_organisms.items():

    pseudo_orgs = sorted(
        org for org in organisms
        if is_pseudomonas(org)
    )

    nonpseudo_orgs = sorted(
        org for org in organisms
        if not is_pseudomonas(org)
    )

    # If an assembly somehow has contradictory taxonomy in IPG,
    # DO NOT automatically delete it.
    if pseudo_orgs and nonpseudo_orgs:

        ambiguous_ipg[assembly] = sorted(organisms)

    elif pseudo_orgs:

        papi29_pseudomonas[assembly] = pseudo_orgs

    elif nonpseudo_orgs:

        papi29_nonpseudomonas[assembly] = nonpseudo_orgs


print()
print("========================================")
print("PAPI29 IPG CLASSIFICATION")
print("========================================")

print(
    f"PAPI29 Pseudomonas assemblies:       "
    f"{len(papi29_pseudomonas)}"
)

print(
    f"PAPI29 non-Pseudomonas assemblies:   "
    f"{len(papi29_nonpseudomonas)}"
)

print(
    f"Ambiguous/conflicting IPG assemblies: "
    f"{len(ambiguous_ipg)}"
)


# ============================================================
# 3. READ CURRENT genomes.txt
# ============================================================

existing_paths = []

with GENOMES_FILE.open() as handle:

    for line in handle:

        line = line.strip()

        if line:
            existing_paths.append(line)


print()
print("Original genomes.txt entries:", len(existing_paths))


# ============================================================
# 4. CLEAN THE CURRENT LIST
# ============================================================
#
# THIS IS THE IMPORTANT RULE:
#
# If an assembly is absent from IPG:
#     KEEP IT.
#
# If IPG explicitly says Pseudomonas:
#     KEEP IT.
#
# If IPG explicitly says non-Pseudomonas:
#     REMOVE IT.
#
# If IPG is contradictory:
#     KEEP IT and report it.
#
# ============================================================

clean_paths = []

removed_nonpseudo = []
retained_background = []
retained_pseudo_papi29 = []
retained_ambiguous = []
unparsed_existing = []


for path in existing_paths:

    assembly = get_accession(path)

    # --------------------------------------------
    # No accession recoverable from filename/path
    # --------------------------------------------

    if assembly is None:

        # Do NOT delete automatically.
        clean_paths.append(path)
        unparsed_existing.append(path)

        continue

    # --------------------------------------------
    # Explicitly NON-Pseudomonas PAPI29 assembly
    # --------------------------------------------

    if assembly in papi29_nonpseudomonas:

        removed_nonpseudo.append(
            (
                assembly,
                "; ".join(papi29_nonpseudomonas[assembly]),
                path
            )
        )

        continue

    # --------------------------------------------
    # Explicitly Pseudomonas PAPI29 assembly
    # --------------------------------------------

    if assembly in papi29_pseudomonas:

        clean_paths.append(path)

        retained_pseudo_papi29.append(
            (
                assembly,
                "; ".join(papi29_pseudomonas[assembly]),
                path
            )
        )

        continue

    # --------------------------------------------
    # Contradictory IPG taxonomy
    # --------------------------------------------

    if assembly in ambiguous_ipg:

        clean_paths.append(path)

        retained_ambiguous.append(
            (
                assembly,
                "; ".join(ambiguous_ipg[assembly]),
                path
            )
        )

        continue

    # --------------------------------------------
    # NOT PRESENT IN IPG
    #
    # This is assumed to belong to the original
    # dereplicated Pseudomonas background.
    #
    # KEEP.
    # --------------------------------------------

    clean_paths.append(path)

    retained_background.append(
        (
            assembly,
            "",
            path
        )
    )


print()
print("========================================")
print("CURRENT GENOME LIST CLEANING")
print("========================================")

print(
    f"Original entries:                       "
    f"{len(existing_paths)}"
)

print(
    f"Explicit non-Pseudomonas PAPI29 removed:"
    f" {len(removed_nonpseudo)}"
)

print(
    f"Pseudomonas PAPI29 entries retained:    "
    f"{len(retained_pseudo_papi29)}"
)

print(
    f"Background genomes absent from IPG kept:"
    f" {len(retained_background)}"
)

print(
    f"Ambiguous IPG entries kept:             "
    f"{len(retained_ambiguous)}"
)

print(
    f"Paths without parseable GCF/GCA kept:   "
    f"{len(unparsed_existing)}"
)


# ============================================================
# 5. INDEX ALL GENOME FASTA FILES UNDER ncbi/
# ============================================================

available_genomes = defaultdict(list)

valid_suffixes = (
    ".fna",
    ".fna.gz",
    ".fa",
    ".fa.gz",
    ".fasta",
    ".fasta.gz",
)


for path in NCBI_DIR.rglob("*"):

    if not path.is_file():
        continue

    path_string = str(path)

    if not path_string.lower().endswith(valid_suffixes):
        continue

    assembly = get_accession(path_string)

    if assembly is not None:
        available_genomes[assembly].append(path)


print()
print(
    "Assemblies with genome FASTA under ncbi/:",
    len(available_genomes)
)


# ============================================================
# 6. FIND ASSEMBLIES ALREADY PRESENT AFTER CLEANING
# ============================================================

present_assemblies = set()

for path in clean_paths:

    assembly = get_accession(path)

    if assembly:
        present_assemblies.add(assembly)


# ============================================================
# 7. ADD MISSING PSEUDOMONAS PAPI29 ASSEMBLIES
# ============================================================

papi29_already_present = []
papi29_added = []
papi29_missing_locally = []


for assembly in sorted(papi29_pseudomonas):

    organisms = "; ".join(
        papi29_pseudomonas[assembly]
    )

    # --------------------------------------------
    # Already represented: DO NOT duplicate
    # --------------------------------------------

    if assembly in present_assemblies:

        papi29_already_present.append(
            (
                assembly,
                organisms,
                ""
            )
        )

        continue

    # --------------------------------------------
    # Not represented but genome exists locally
    # --------------------------------------------

    if assembly in available_genomes:

        candidates = available_genomes[assembly]

        # Prefer files containing "genomic.fna" if NCBI has
        # left multiple FASTA files for the assembly.
        preferred = [
            p for p in candidates
            if "genomic.fna" in p.name
        ]

        if preferred:
            genome_path = sorted(preferred)[0]
        else:
            genome_path = sorted(candidates)[0]

        clean_paths.append(str(genome_path))

        present_assemblies.add(assembly)

        papi29_added.append(
            (
                assembly,
                organisms,
                str(genome_path)
            )
        )

    # --------------------------------------------
    # Pseudomonas PAPI29 assembly exists in IPG
    # but we don't have its genome locally.
    # --------------------------------------------

    else:

        papi29_missing_locally.append(
            (
                assembly,
                organisms,
                ""
            )
        )


# ============================================================
# 8. FINAL DEDUPLICATION
# ============================================================
#
# Deduplicate by assembly accession, NOT filename.
#
# ============================================================

final_paths = []

seen_assemblies = set()
seen_unparsed = set()

duplicates_removed = []


for path in clean_paths:

    assembly = get_accession(path)

    # Genome with identifiable GCF/GCA
    if assembly:

        if assembly in seen_assemblies:

            duplicates_removed.append(
                (
                    assembly,
                    "",
                    path
                )
            )

            continue

        seen_assemblies.add(assembly)
        final_paths.append(path)

    # Weird/unparseable filename: deduplicate only exact path
    else:

        if path in seen_unparsed:
            continue

        seen_unparsed.add(path)
        final_paths.append(path)


# ============================================================
# 9. WRITE FINAL GENOME LIST
# ============================================================

with OUTPUT.open("w") as handle:

    for path in final_paths:
        handle.write(path + "\n")


# ============================================================
# 10. WRITE FULL AUDIT REPORT
# ============================================================

with REPORT.open("w") as handle:

    handle.write(
        "status\tassembly\torganism\tpath\n"
    )

    # Things actually removed
    for assembly, organism, path in removed_nonpseudo:

        handle.write(
            f"REMOVED_NON_PSEUDOMONAS_PAPI29\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # PAPI29+ Pseudomonas already in original list
    for assembly, organism, path in retained_pseudo_papi29:

        handle.write(
            f"RETAINED_PSEUDOMONAS_PAPI29\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # Background genomes not represented in IPG
    for assembly, organism, path in retained_background:

        handle.write(
            f"RETAINED_BACKGROUND_NOT_IN_IPG\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # Ambiguous IPG taxonomy
    for assembly, organism, path in retained_ambiguous:

        handle.write(
            f"RETAINED_AMBIGUOUS_IPG_TAXONOMY\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # Pseudomonas PAPI29 already represented after cleaning
    for assembly, organism, path in papi29_already_present:

        handle.write(
            f"PSEUDOMONAS_PAPI29_ALREADY_PRESENT\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # Newly added Pseudomonas PAPI29 genomes
    for assembly, organism, path in papi29_added:

        handle.write(
            f"ADDED_PSEUDOMONAS_PAPI29\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # Pseudomonas PAPI29 hits whose genomes are unavailable locally
    for assembly, organism, path in papi29_missing_locally:

        handle.write(
            f"PSEUDOMONAS_PAPI29_NOT_FOUND_LOCALLY\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # Duplicate assembly paths
    for assembly, organism, path in duplicates_removed:

        handle.write(
            f"DUPLICATE_ASSEMBLY_REMOVED\t"
            f"{assembly}\t"
            f"{organism}\t"
            f"{path}\n"
        )

    # Paths for which accession couldn't be parsed
    for path in unparsed_existing:

        handle.write(
            f"RETAINED_UNPARSED_PATH\t"
            f"\t\t"
            f"{path}\n"
        )


# ============================================================
# 11. FINAL SANITY CHECK
# ============================================================

# There should be ZERO assemblies in the final list that are
# explicitly classified as non-Pseudomonas by the IPG.

bad_final = []

for path in final_paths:

    assembly = get_accession(path)

    if assembly in papi29_nonpseudomonas:

        bad_final.append(
            (
                assembly,
                "; ".join(
                    papi29_nonpseudomonas[assembly]
                ),
                path
            )
        )


# ============================================================
# 12. SUMMARY
# ============================================================

print()
print("========================================")
print("FINAL SUMMARY")
print("========================================")

print(
    f"Original genome-list entries:             "
    f"{len(existing_paths)}"
)

print(
    f"Explicit non-Pseudomonas PAPI29 removed:   "
    f"{len(removed_nonpseudo)}"
)

print(
    f"Pseudomonas PAPI29 already represented:    "
    f"{len(papi29_already_present)}"
)

print(
    f"Pseudomonas PAPI29 newly added:             "
    f"{len(papi29_added)}"
)

print(
    f"Pseudomonas PAPI29 missing locally:         "
    f"{len(papi29_missing_locally)}"
)

print(
    f"Duplicate assembly entries removed:         "
    f"{len(duplicates_removed)}"
)

print(
    f"Final genome-list entries:                  "
    f"{len(final_paths)}"
)

print(
    f"Known non-Pseudomonas IPG assemblies left:  "
    f"{len(bad_final)}"
)

print()


if bad_final:

    print("WARNING!")
    print(
        "Known non-Pseudomonas PAPI29 assemblies "
        "remain in the final list:"
    )

    for assembly, organism, path in bad_final:
        print(
            f"  {assembly}\t{organism}\t{path}"
        )

else:

    print(
        "PASS: no assembly explicitly classified as "
        "non-Pseudomonas in the IPG remains in the final list."
    )


print()
print(f"Final GToTree genome list: {OUTPUT}")
print(f"Full audit report:         {REPORT}")
