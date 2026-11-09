# Culsans bioinformatics analysis

This directory contains the scripts, inputs, intermediate files, and reproducible outputs used for the Culsans study. Each analysis script is documented below as the workflow develops.

## Project layout

```text
PAPI29-AllBoinfo/
├── README.md
├── blast_hits.tsv
├── etnav3.py
├── flankguard.py
├── results/                 # etnav3.py output
└── flankguard_results/      # flankguard.py output
```

- `blast_hits.tsv`: BLAST hit table used as input.
- `etnav3.py`: retrieves RefSeq IPG mappings and NCBI source metadata.
- `flankguard.py`: profiles defence and anti-defence systems in representative genomic neighbourhoods.

## Script 1: `etnav3.py`

### Purpose

`etnav3.py` converts a BLAST protein-hit table into a nucleotide-level RefSeq dataset enriched with NCBI Nuccore and BioSample metadata. It helps describe the taxonomic distribution and reported isolation context of sequences matching the query protein.

The script:

1. Reads a BLAST table delimited by tabs, commas, semicolons, or pipes.
2. Extracts protein accessions from column 2 and deduplicates them in input order.
3. Queries NCBI Identical Protein Groups (IPG) for each protein.
4. Retains RefSeq nucleotide mappings and available CDS, organism, strain, assembly, and protein metadata.
5. Retrieves corresponding Nuccore and BioSample source metadata.
6. Derives genus from the reported organism name.
7. Applies keyword-assisted screening labels for clinical status and isolation-source category.
8. Writes TSV/JSON results and an HTML report with an SVG summary figure.

### Requirements and use

`etnav3.py` requires Python 3.9 or newer, internet access unless the cache is complete, and a contact email for NCBI E-utilities. It uses only the Python standard library. An optional NCBI API key can be passed with `--api-key` or `NCBI_API_KEY`.

Run from the project directory, replacing the example email:

```bash
python3 etnav3.py blast_hits.tsv results \
  --email your.email@example.org \
  --header yes
```

`--header yes` is required for the current table because its first fields are `input_line` and `protein_accession`, a pair not reliably recognized by this script's automatic header detection.

Small test run:

```bash
python3 etnav3.py blast_hits.tsv test_results \
  --email your.email@example.org \
  --header yes \
  --max-proteins 5 \
  --max-nuccore 10
```

Cached/offline rerun:

```bash
python3 etnav3.py blast_hits.tsv results \
  --email your.email@example.org \
  --header yes \
  --offline
```

### Main outputs

| Output | Description |
|---|---|
| `blast_hits.tsv` | Normalized selected BLAST hits. |
| `ipg_refseq_hits.tsv` | Protein-to-RefSeq nucleotide mappings returned by IPG. |
| `nuccore_metadata.tsv` | One row per nucleotide accession with Nuccore/BioSample metadata. |
| `refseq_hits_enriched.tsv` | IPG mappings joined to nucleotide metadata. |
| `genus_counts.tsv` | Counts and percentages by inferred genus. |
| `clinical_status_counts.tsv` | Counts and percentages by clinical-status screening label. |
| `source_category_counts.tsv` | Counts and percentages by isolation-source screening category. |
| `failures.tsv` | Failed accessions, stages, and error messages. |
| `summary.json` | Run-level counts and provenance. |
| `hit_summary.svg` | Summary plots. |
| `report.html` | Human-readable summary report. |
| `cache/` | Resumable raw NCBI XML responses. |

Counts use unique RefSeq nucleotide accessions. Automated clinical/source categories are screening aids rather than validated clinical annotations. Review the evidence and original NCBI records before biological interpretation.

## Script 2: FlankGuard (`flankguard.py`)

### Purpose

FlankGuard is a reusable genomic-neighbourhood defence profiler. For each unique protein accession in a BLAST table, it selects one RefSeq nucleotide representative, extracts the target CDS plus 30 kb upstream and 30 kb downstream, runs DefenseFinder with both defence and AntiDefenseFinder models, runs PADLOC, combines the calls, and summarizes how often the query protein occurs in a neighbourhood containing a detected defence or anti-defence system.

For this study, the query proteins are Culsans hits. The program itself does not hard-code the Culsans name and can therefore be reused for other proteins.

### Workflow

1. Detects the protein-accession column (`protein_accession`, `sseqid`, related names, or a user-specified column).
2. Deduplicates RefSeq protein accessions while retaining their input order.
3. Retrieves each protein's NCBI IPG report.
4. Selects one deterministic RefSeq nucleotide representative per protein.
5. Retrieves a genomic window containing the target CDS and the requested flank on both sides.
6. Creates nucleotide FASTA, ordered protein FASTA, GFF3, and gene-map files for every neighbourhood.
7. Runs DefenseFinder with `--antidefensefinder` so defence and anti-defence models are both searched.
8. Runs PADLOC on the same ordered proteins and GFF coordinates.
9. Combines system calls by representative neighbourhood.
10. Writes association counts, system frequencies, a frequency plot, and a presence/absence heatmap.

### Representative-selection rule

Selection is deterministic and recorded in `representatives.tsv`. FlankGuard first requires usable CDS coordinates, then prefers nucleotide accessions in this order:

1. curated RefSeq records beginning `NC_`;
2. RefSeq complete-replicon records beginning `NZ_CP`;
3. `NZ_AP` or `NZ_BA` records;
4. other `NZ_` records;
5. other RefSeq mappings.

Within a tier, mappings with assembly, organism, and strain metadata are preferred, followed by a lexical accession tie-break. This produces one reproducible context per protein but does not claim that the selected record is biologically representative of every genome containing that identical protein.

### Replicon edges

FlankGuard retrieves the replicon length and topology from NCBI. Windows on circular replicons wrap across the coordinate origin when necessary. Windows on linear records are truncated at the sequence boundary and explicitly marked in the `edge_status` column. The exact retrieved interval or intervals are recorded in `neighborhood_manifest.tsv`.

### Existing environments

DefenseFinder and PADLOC are already installed locally in Conda environments named:

- `defensefinder`
- `padloc`

These are the script defaults. Alternative environment names can be supplied with `--defensefinder-env` and `--padloc-env`. FlankGuard itself uses only the Python standard library.

Environment validation on 13 September 2026 confirmed that DefenseFinder 2.0.1 runs successfully with defence and anti-defence models. The `padloc` environment and its R dependencies are also installed. FlankGuard places each selected Conda environment at the front of the executable path so PADLOC uses the environment's matching R installation even when the workflow is launched from a desktop application or IDE. It checks both scanner environments before making NCBI requests, so missing dependencies are reported early.

### Basic use

Replace the example email address with the address to send to NCBI:

```bash
python3 flankguard.py blast_hits.tsv flankguard_results \
  --email your.email@example.org \
  --defensefinder-env defensefinder \
  --padloc-env padloc \
  --jobs 4 \
  --threads 1 \
  --heatmap-color '#264979'
```

The current `blast_hits.tsv` contains a `protein_accession` field, which FlankGuard detects automatically. A different field can be selected explicitly:

```bash
python3 flankguard.py hits.tsv flankguard_results \
  --email your.email@example.org \
  --protein-column sseqid
```

### Test before the complete run

A small run is strongly recommended before processing the complete hit table:

```bash
python3 flankguard.py blast_hits.tsv flankguard_test \
  --email your.email@example.org \
  --max-proteins 3 \
  --threads 2
```

`--dry-run` retrieves and prepares contexts and prints the scanner commands without executing DefenseFinder or PADLOC. `--offline` makes no NCBI requests and succeeds only when all required XML files already exist in the selected output directory's cache.

### Important arguments

| Argument | Meaning |
|---|---|
| `--email EMAIL` | Contact email required by NCBI E-utilities. |
| `--api-key KEY` | Optional NCBI API key; can also be supplied as `NCBI_API_KEY`. |
| `--protein-column NAME` | Explicit protein-accession column; otherwise detected automatically. |
| `--flank-bp N` | Bases retrieved on each side of the target CDS. Default: `30000`. |
| `--defensefinder-env NAME` | Conda environment containing DefenseFinder. Default: `defensefinder`. |
| `--padloc-env NAME` | Conda environment containing PADLOC. Default: `padloc`. |
| `--conda PATH` | Conda/Mamba executable. Automatically detected when possible. |
| `--threads N` | CPUs supplied to each scanner. Default: `1`. |
| `--jobs N` | Independent neighbourhoods scanned concurrently. Default: `1`. Total potential CPU use is approximately `jobs × threads`. |
| `--heatmap-color HEX` | Main heatmap and plot colour. Default: `#264979`. |
| `--heatmap-top-systems N` | Maximum number of frequent system columns plotted. Default: `25`. |
| `--heatmap-max-neighborhoods N` | Maximum number of neighbourhood rows plotted. Default: `250`. |
| `--max-proteins N` | Processes only the first `N` unique proteins for testing. |
| `--dry-run` | Prepares inputs and commands but does not execute the scanners. |
| `--skip-tool-checks` | Skips scanner preflight checks; use only when checks are unavailable but the environments are known to work. |
| `--offline` | Uses cached NCBI XML only. |
| `--skip-defensefinder` | Skips DefenseFinder execution and parsing. |
| `--skip-padloc` | Skips PADLOC execution and parsing. |

Quote hexadecimal colours on the command line so the shell does not interpret `#` as the start of a comment.

### Outputs

| Output | Description |
|---|---|
| `representatives.tsv` | Chosen RefSeq mapping, candidate count, and selection rule for each protein. |
| `neighborhood_manifest.tsv` | Exact nucleotide context, coordinates, topology, length, edge handling, and quality-control fields. |
| `neighborhoods/<protein>/neighborhood.fna` | Extracted nucleotide window. |
| `neighborhoods/<protein>/neighborhood.faa` | Ordered translated CDS features used by both scanners. |
| `neighborhoods/<protein>/neighborhood.gff` | Matching GFF3 coordinates used by PADLOC. |
| `neighborhoods/<protein>/genes.tsv` | Map between local feature IDs and original NCBI annotations. |
| `neighborhoods/<protein>/defensefinder/` | Raw DefenseFinder and AntiDefenseFinder results plus a run log. |
| `neighborhoods/<protein>/padloc/` | Raw PADLOC results plus a run log. |
| `combined_systems.tsv` | Standardized combined system calls with the original result retained as JSON. |
| `neighborhood_summary.tsv` | One row per protein representative, detected systems, and defence-environment status. |
| `query_environment_counts.tsv` | Number and percentage of query-protein neighbourhoods with or without detected systems. For this dataset, these are Culsans neighbourhood counts. |
| `system_frequencies.tsv` | Number and percentage of representative neighbourhoods containing each original tool-specific system call (provenance table). |
| `consensus_system_frequencies.tsv` | Frequencies of harmonised, nonredundant biological systems after overlapping cross-program calls are merged. |
| `system_frequencies.svg` | Main horizontal bar chart using harmonised, nonredundant system loci. |
| `tool_specific_system_frequencies.svg` | Provenance plot retaining DefenseFinder and PADLOC calls separately. |
| `defence_environment_percentage.svg` | Percentage plot comparing Culsans neighbourhoods with and without detected systems. |
| `system_heatmap.svg` | Multi-column neighbourhood-by-system presence/absence heatmap with square outlined cells. |
| `system_heatmap_A4.svg` | A4-landscape genus-by-system prevalence heatmap for a manuscript figure. |
| `association_breadth.svg` | Compact summary of the prevalence and taxonomic breadth of defence association. |
| `system_presence_matrix.tsv` | Reusable 0/1 neighbourhood-by-system matrix for plotting a heatmap in R, Python, or spreadsheet software. |
| `tool_specific_system_frequency_strip.svg` | Optional provenance-only frequency strip with calls separated by program. |
| `supplementary_system_calls.tsv` | One row per DefenseFinder, AntiDefenseFinder, or PADLOC system call with constituent genes, RefSeq protein accessions, coordinates, strands, distance, and position relative to Culsans. |
| `supplementary_consensus_loci.tsv` | One row per nonredundant biological locus, retaining both programs' original labels and all supporting genes and coordinates. |
| `supplementary_neighborhoods.tsv` | One combined row per Culsans neighbourhood with separate and integrated scanner results. |
| `Culsans_defence_supplement.xlsx` | Formatted manuscript supplementary workbook containing the combined neighbourhood and full system-call tables. |
| `failures.tsv` | Failures from IPG retrieval, context extraction, DefenseFinder, or PADLOC. |
| `summary.json` | Run parameters and headline counts. |
| `cache/` | Raw cached IPG, Nuccore summary, and region XML responses. |

The main counting unit is one selected RefSeq nucleotide neighbourhood per unique protein accession. Original tool-prefixed calls are retained for provenance. A second consensus layer merges calls at the same genomic locus when they share constituent genes, or when their genomic spans overlap and harmonised system names agree. Naming variants are normalised before comparison—for example, `CBASS_I` and `cbass_type_I` become CBASS type I, and `RM_Type_II` and `RM_type_II` become RM type II. The original program names remain visible in the detailed tables.

`DMS_other` is excluded from the biological consensus, headline association, plots, and merged plotting tables. It remains only in the raw tool-specific call table so the original PADLOC output is auditable.

The association percentage is calculated as the number of analysable Culsans representative neighbourhoods containing at least one retained DefenseFinder, AntiDefenseFinder, or PADLOC call divided by all analysable representative neighbourhoods. `DMS_other` is excluded from the merged biological analysis; it remains only in the raw program-call table for auditability. The 17 proteins without a RefSeq IPG mapping are reported as exclusions and are not included in this denominator.



