# Culsans--LPS association workflow

This workflow tests whether LPS, lipid A, core-oligosaccharide and
O-antigen/nucleotide-sugar KOs are associated with Culsans distribution
across a *Pseudomonas* genome cohort, while accounting for population
structure.

## Workflow overview

The final analysis was performed in the following order:

``` text
Original Pseudomonas genome set + Culsans IPG distribution (Culsansdistri)
                              |
                              v
                    fix_Culsans_genomes.py
                              |
             +----------------+----------------+
             |                                 |
             v                                 v
  genomes_merged_gotree.txt          cleaned Pseudomonas cohort
             |                                 |
             v                                 v
          GToTree                    faa_pseudomonas.txt
             |                                 |
             v                                 v
      Aligned_SCGs.faa                   KoFamScan.py
             |                                 |
             v                                 v
    population phylogeny            ko_presence_unambiguous.tsv
             |                                 |
             +---------------+-----------------+
                             |
                             v
                 lps-Culsans-association.py
                             |
                             v
             phylogeny-adjusted KO associations
```



## 1. Construct the Pseudomonas/Culsans genome cohort --- `fix_Culsans_genomes.py`

The starting genome set contained *Pseudomonas* genomes plus
representative genomes added from the Culsans IPG distribution. Because
the downstream biological work concerns Culsans in *Pseudomonas*, Culsans
representatives belonging to non-*Pseudomonas* genera were removed.
Missing *Pseudomonas* Culsans representatives were added while avoiding
duplicate assembly accessions.

The resulting cleaned/merged cohort was used consistently for
phylogenetic reconstruction and LPS/pili annotation.

Important downstream cohort files include:

-   `genomes_merged_gotree.txt` --- genome inputs used for GToTree.
-   `gtotree_merged_accessions.txt` / the subsequently cleaned
    Pseudomonas accession list --- assembly accessions representing the
    analysis cohort.
-   `faa_pseudomonas.txt` --- paths to protein FASTA files for the same
    Pseudomonas cohort, used by `KoFamScan.py`.

## 2. Population phylogeny --- GToTree

A phylogenomic alignment was generated from the merged genome cohort
using GToTree and the Gammaproteobacteria single-copy-gene marker set:

``` bash
GToTree \
  -f genomes_merged_gotree.txt \
  -H Gammaproteobacteria \
  -j 8 \
  -o GToTree_Culsans
```


## 3. Targeted LPS KOfam annotation --- `KoFamScan.py`

Protein FASTA files corresponding to the cleaned Pseudomonas cohort were
collected in:

`faa_pseudomonas.txt`

`KoFamScan.py` runs targeted KofamScan against the LPS/O-antigen KO
sets using calibrated KOfam thresholds. Existing gene/product
annotations are not used to assign the targeted KOs.

Example:

``` bash
python KoFamScan.py \
  --root . \
  --faa-list faa_pseudomonas.txt \
  --cpu 8 \
  --outdir lps_kofam_new
```

The principal downstream matrix is:

`lps_kofam_new/ko_presence_unambiguous.tsv`

A value of `1` means that at least one unambiguous protein hit passed
the KO-specific KOfam threshold. A value of `0` means that no qualifying
hit was detected; it should not be interpreted as proof of biological
absence.

## 4. Culsans--LPS association --- `lps-Culsans-association.py`

The final script combines four components:

1.  the cleaned *Pseudomonas* assembly cohort;
2.  Culsans distribution from `Culsansdistri`;
3.  LPS/O-antigen KO presence from `ko_presence_unambiguous.tsv`; and
4.  the population phylogeny derived from the GToTree SCG alignment.

Example:

``` bash
python lps-Culsans-association.py \
  --tree <final_population_tree> \
  --kofam lps_kofam_new/ko_presence_unambiguous.tsv \
  --ipg Culsansdistri \
  --cohort gtotree_merged_accessions_pseudomonas.txt \
  --outdir Culsans_LPS_association \
  --palette '#F5EBDD' '#F2765E' '#315B8C' '#413333' \
  --phylo-variance 0.80 \
  --top-n 20
```

Culsans status is assigned from the supplied IPG table: cohort
*Pseudomonas* assemblies represented in the Culsans IPG are Culsans+,
while other cohort assemblies are Culsans−.

For each KO, the script reports Culsans+/− prevalence, Fisher's exact
test, odds ratio and Benjamini--Hochberg FDR. Population structure is
represented by principal coordinates calculated from pairwise patristic
distances in the supplied tree. Enough coordinates are retained to meet
`--phylo-variance`. Each KO is then tested by comparing logistic models
containing the phylogenetic coordinates alone with models additionally
containing KO presence using a likelihood-ratio test. Phylogeny-adjusted
P values are Benjamini--Hochberg corrected.

## Main association outputs

-   `lps_association_results.tsv` --- complete results for all tested
    KOs.
-   `significant_phylogeny_adjusted_associations.tsv` --- KOs with
    phylogeny-adjusted FDR \< 0.05.
-   `naive_only_associations.tsv` --- Fisher-significant KOs that do not
    remain significant after phylogenetic adjustment.
-   `phylogenetic_coordinates.tsv` and `phylogenetic_eigenvalues.tsv`
    --- population-structure covariates and represented variance.
-   `01_association_scatter.*` --- association effect versus
    phylogeny-adjusted significance.
-   `02_prevalence_heatmap.*` --- prevalence of selected KOs in Culsans−
    and Culsans+ genomes.
-   `03_tree_KO_heatmap.*` --- phylogenetic distribution of Culsans and
    selected KOs.

## Interpretation

For phylogeny-adjusted Culsans enrichment, use **`phylo_FDR < 0.05` and
`phylo_beta > 0`**. For depletion, use **`phylo_FDR < 0.05` and
`phylo_beta < 0`**.

Raw prevalence direction and the phylogenetically adjusted coefficient
can occasionally disagree because the adjusted model estimates the KO
association conditional on population structure.

Very common KOs can produce large odds ratios from a small number of
KO-negative genomes. Therefore, prevalence difference, adjusted effect
size, pathway-level coherence and biological function should be
considered alongside P values and FDR when prioritising candidate LPS
pathways for experimental testing.
