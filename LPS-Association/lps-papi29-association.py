#!/usr/bin/env python3

"""
PAPI29_LPS_Association.py

Association analysis between PAPI29/Culsans presence and LPS/O-antigen
KOfam annotations while accounting for phylogenetic population structure.

PAPI29 status:
    Pseudomonas assemblies present in the supplied IPG file are PAPI29+.
    Cohort assemblies absent from the IPG are PAPI29-.

Analyses:
    1. Fisher's exact test for every KO.
    2. Odds ratio and prevalence difference.
    3. Benjamini-Hochberg correction of naive P values.
    4. Principal-coordinate representation of patristic tree distances.
    5. Logistic regression controlling for phylogenetic coordinates.
    6. Likelihood-ratio test:
           PAPI29 ~ phylogeny
       versus
           PAPI29 ~ phylogeny + KO
    7. Benjamini-Hochberg correction of phylogeny-adjusted P values.

Figures:
    1. Effect-size/significance scatter.
    2. PAPI29+ versus PAPI29- KO prevalence heatmap.
    3. Phylogenetic tree with PAPI29 status and selected KO presence.

Example:

python PAPI29_LPS_Association.py \
    --tree GToTree_PAPI29/GToTree_PAPI29.tre \
    --kofam lps_kofam_new/ko_presence_unambiguous.tsv \
    --ipg papi29distri \
    --cohort gtotree_merged_accessions.txt \
    --outdir PAPI29_LPS_association \
    --palette '#F5EBDD' '#F2765E' '#315B8C' '#413333' \
    --phylo-variance 0.80 \
    --top-n 20
"""

import argparse
import csv
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from scipy.stats import fisher_exact, chi2

import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning

from Bio import Phylo

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


warnings.filterwarnings(
    "ignore",
    category=PerfectSeparationWarning
)


# ============================================================
# ACCESSION PARSING
# ============================================================

ACC_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def accession(text):

    m = ACC_RE.search(str(text))

    if m:
        return m.group(1)

    return None


# ============================================================
# BH FDR
# ============================================================

def bh_fdr(values):

    values = np.asarray(values, dtype=float)

    result = np.full(
        len(values),
        np.nan
    )

    valid = np.isfinite(values)

    p = values[valid]

    if len(p) == 0:
        return result

    order = np.argsort(p)
    ranked = p[order]

    n = len(ranked)

    adjusted = (
        ranked
        * n
        / np.arange(1, n + 1)
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.minimum(
        adjusted,
        1.0
    )

    tmp = np.empty(n)

    tmp[order] = adjusted

    result[valid] = tmp

    return result


# ============================================================
# COHORT
# ============================================================

def read_cohort(path):

    assemblies = []

    seen = set()

    for line in Path(path).read_text().splitlines():

        line = line.strip()

        if not line:
            continue

        acc = accession(line)

        if acc is None:
            raise ValueError(
                "Could not extract GCF/GCA accession "
                f"from cohort line:\n{line}"
            )

        if acc in seen:
            raise ValueError(
                f"Duplicate cohort accession: {acc}"
            )

        seen.add(acc)
        assemblies.append(acc)

    if not assemblies:
        raise ValueError(
            "Cohort file contains no assemblies."
        )

    return assemblies


# ============================================================
# PAPI29 STATUS FROM IPG
# ============================================================

def read_papi29_ipg(path):

    positive = set()

    with Path(path).open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        required = {
            "assembly",
            "organism"
        }

        if not required.issubset(
            reader.fieldnames or []
        ):
            raise ValueError(
                "IPG file must contain columns "
                "'assembly' and 'organism'."
            )

        for row in reader:

            org = (
                row["organism"]
                .strip()
                .lower()
            )

            acc = accession(
                row["assembly"]
            )

            if (
                acc
                and org.startswith(
                    "pseudomonas "
                )
            ):
                positive.add(acc)

    return positive


# ============================================================
# KOFAM MATRIX
# ============================================================

def read_kofam(path):

    df = pd.read_csv(
        path,
        sep="\t"
    )

    if "assembly_accession" not in df.columns:
        raise ValueError(
            "KOfam matrix has no "
            "'assembly_accession' column."
        )

    df["assembly_accession"] = (
        df["assembly_accession"]
        .astype(str)
        .map(accession)
    )

    if df["assembly_accession"].isna().any():

        bad = df[
            df["assembly_accession"].isna()
        ]

        raise ValueError(
            "Could not parse some KOfam "
            "assembly accessions:\n"
            + str(bad.head())
        )

    if df[
        "assembly_accession"
    ].duplicated().any():

        duplicates = df.loc[
            df[
                "assembly_accession"
            ].duplicated(),
            "assembly_accession"
        ].tolist()

        raise ValueError(
            "Duplicate assemblies in KOfam "
            f"matrix: {duplicates[:10]}"
        )

    # This column was retained by the KOfam
    # annotation script but isn't used here.
    if "culsans_positive" in df.columns:

        df = df.drop(
            columns=["culsans_positive"]
        )

    ko_columns = [
        c for c in df.columns
        if c != "assembly_accession"
    ]

    if not ko_columns:
        raise ValueError(
            "No KO columns found."
        )

    for ko in ko_columns:

        df[ko] = pd.to_numeric(
            df[ko],
            errors="raise"
        )

        bad_values = set(
            df[ko].dropna().unique()
        ) - {0, 1}

        if bad_values:
            raise ValueError(
                f"{ko} contains values other "
                f"than 0/1: {bad_values}"
            )

    return df, ko_columns


# ============================================================
# TREE TERMINALS
# ============================================================

def get_tree_terminals(tree):

    terminals = {}

    for terminal in tree.get_terminals():

        acc = accession(
            terminal.name
        )

        if acc:

            if acc in terminals:
                raise ValueError(
                    "Duplicate assembly in tree: "
                    + acc
                )

            terminals[acc] = terminal

    return terminals


# ============================================================
# PATRISTIC DISTANCE MATRIX
# ============================================================

def tree_distance_matrix(
    tree,
    assemblies,
    terminals
):

    n = len(assemblies)

    D = np.zeros(
        (n, n),
        dtype=float
    )

    print(
        "\nCalculating patristic distance matrix...",
        flush=True
    )

    for i in range(n):

        if (
            i == 0
            or (i + 1) % 100 == 0
            or i == n - 1
        ):
            print(
                f"  {i + 1}/{n}",
                flush=True
            )

        ti = terminals[
            assemblies[i]
        ]

        for j in range(
            i + 1,
            n
        ):

            tj = terminals[
                assemblies[j]
            ]

            d = tree.distance(
                ti,
                tj
            )

            D[i, j] = d
            D[j, i] = d

    return D


# ============================================================
# PHYLOGENETIC PRINCIPAL COORDINATES
# ============================================================

def phylogenetic_coordinates(
    D,
    variance_target
):

    n = D.shape[0]

    J = (
        np.eye(n)
        - np.ones(
            (n, n)
        ) / n
    )

    B = (
        -0.5
        * J
        @ (D ** 2)
        @ J
    )

    eigenvalues, eigenvectors = (
        np.linalg.eigh(B)
    )

    order = np.argsort(
        eigenvalues
    )[::-1]

    eigenvalues = (
        eigenvalues[order]
    )

    eigenvectors = (
        eigenvectors[:, order]
    )

    positive = (
        eigenvalues > 1e-10
    )

    eigenvalues = (
        eigenvalues[positive]
    )

    eigenvectors = (
        eigenvectors[:, positive]
    )

    if len(eigenvalues) == 0:
        raise ValueError(
            "Tree produced no positive "
            "phylogenetic eigenvalues."
        )

    explained = (
        eigenvalues
        / eigenvalues.sum()
    )

    cumulative = np.cumsum(
        explained
    )

    n_axes = (
        np.searchsorted(
            cumulative,
            variance_target
        )
        + 1
    )

    coords = (
        eigenvectors[:, :n_axes]
        * np.sqrt(
            eigenvalues[:n_axes]
        )
    )

    print(
        "\nPhylogenetic population structure:"
    )

    print(
        "  Positive eigenvalues:",
        len(eigenvalues)
    )

    print(
        "  Axes retained:",
        n_axes
    )

    print(
        "  Variance represented:",
        f"{cumulative[n_axes - 1]:.3f}"
    )

    return (
        coords,
        eigenvalues,
        explained,
        n_axes
    )


# ============================================================
# LOGISTIC LIKELIHOOD-RATIO TEST
# ============================================================

def phylogenetic_logistic_lrt(
    y,
    ko,
    phylo
):

    # Null:
    # PAPI29 ~ phylogenetic coordinates

    X0 = sm.add_constant(
        phylo,
        has_constant="add"
    )

    # Full:
    # PAPI29 ~ phylogenetic coordinates + KO

    X1 = np.column_stack([
        phylo,
        ko
    ])

    X1 = sm.add_constant(
        X1,
        has_constant="add"
    )

    try:

        null = sm.Logit(
            y,
            X0
        ).fit(
            disp=False,
            maxiter=200
        )

        full = sm.Logit(
            y,
            X1
        ).fit(
            disp=False,
            maxiter=200
        )

        LR = max(
            0.0,
            2.0 * (
                full.llf
                - null.llf
            )
        )

        p = chi2.sf(
            LR,
            df=1
        )

        beta = (
            full.params[-1]
        )

        adjusted_or = (
            np.exp(beta)
            if np.isfinite(beta)
            else np.nan
        )

        return (
            beta,
            adjusted_or,
            LR,
            p,
            True
        )

    except Exception:

        return (
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            False
        )


# ============================================================
# ASSOCIATION TESTS
# ============================================================

def run_associations(
    matrix,
    ko_columns,
    phylo
):

    y = matrix[
        "PAPI29"
    ].to_numpy(
        dtype=int
    )

    n_pos = int(
        y.sum()
    )

    n_neg = int(
        len(y) - n_pos
    )

    results = []

    print(
        "\nTesting KO associations...",
        flush=True
    )

    for number, ko_name in enumerate(
        ko_columns,
        start=1
    ):

        if (
            number == 1
            or number % 25 == 0
            or number == len(ko_columns)
        ):
            print(
                f"  {number}/{len(ko_columns)}",
                flush=True
            )

        x = matrix[
            ko_name
        ].to_numpy(
            dtype=int
        )

        # ----------------------------------------
        # 2 x 2 table
        #
        #             KO+   KO-
        # PAPI29+      a     b
        # PAPI29-      c     d
        # ----------------------------------------

        a = int(
            ((y == 1) & (x == 1)).sum()
        )

        b = int(
            ((y == 1) & (x == 0)).sum()
        )

        c = int(
            ((y == 0) & (x == 1)).sum()
        )

        d = int(
            ((y == 0) & (x == 0)).sum()
        )

        prevalence_pos = (
            a / n_pos
            if n_pos
            else np.nan
        )

        prevalence_neg = (
            c / n_neg
            if n_neg
            else np.nan
        )

        prevalence_difference = (
            prevalence_pos
            - prevalence_neg
        )

        # Fisher
        fisher_or, fisher_p = (
            fisher_exact(
                [[a, b], [c, d]],
                alternative="two-sided"
            )
        )

        # Haldane-Anscombe corrected OR,
        # useful for plotting finite log2 ORs
        corrected_or = (
            ((a + 0.5) * (d + 0.5))
            /
            ((b + 0.5) * (c + 0.5))
        )

        log2_or = np.log2(
            corrected_or
        )

        # ----------------------------------------
        # Phylogeny-adjusted model
        # ----------------------------------------

        (
            beta,
            phylo_or,
            lr,
            phylo_p,
            model_ok
        ) = phylogenetic_logistic_lrt(
            y,
            x,
            phylo
        )

        results.append({
            "KO":
                ko_name,

            "PAPI29_pos_KO_pos":
                a,

            "PAPI29_pos_KO_neg":
                b,

            "PAPI29_neg_KO_pos":
                c,

            "PAPI29_neg_KO_neg":
                d,

            "prevalence_PAPI29_pos":
                prevalence_pos,

            "prevalence_PAPI29_neg":
                prevalence_neg,

            "prevalence_difference":
                prevalence_difference,

            "fisher_odds_ratio":
                fisher_or,

            "haldane_corrected_OR":
                corrected_or,

            "log2_OR":
                log2_or,

            "fisher_p":
                fisher_p,

            "phylo_beta":
                beta,

            "phylo_odds_ratio":
                phylo_or,

            "phylo_LR":
                lr,

            "phylo_p":
                phylo_p,

            "phylo_model_converged":
                model_ok
        })

    results = pd.DataFrame(
        results
    )

    results[
        "fisher_FDR"
    ] = bh_fdr(
        results[
            "fisher_p"
        ].to_numpy()
    )

    results[
        "phylo_FDR"
    ] = bh_fdr(
        results[
            "phylo_p"
        ].to_numpy()
    )

    results[
        "naive_significant"
    ] = (
        results[
            "fisher_FDR"
        ] < 0.05
    )

    results[
        "phylo_significant"
    ] = (
        results[
            "phylo_FDR"
        ] < 0.05
    )

    results[
        "association_direction"
    ] = np.where(
        results[
            "prevalence_difference"
        ] > 0,
        "enriched_with_PAPI29",
        np.where(
            results[
                "prevalence_difference"
            ] < 0,
            "depleted_with_PAPI29",
            "no_difference"
        )
    )

    return results


# ============================================================
# FIGURE SAVING
# ============================================================

def save_figure(
    fig,
    outdir,
    basename
):

    fig.savefig(
        outdir / (
            basename + ".pdf"
        ),
        bbox_inches="tight"
    )

    fig.savefig(
        outdir / (
            basename + ".svg"
        ),
        bbox_inches="tight"
    )

    fig.savefig(
        outdir / (
            basename + ".png"
        ),
        dpi=400,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# FIGURE 1:
# EFFECT SIZE / SIGNIFICANCE
# ============================================================

def association_scatter(
    results,
    outdir,
    palette
):

    light, positive, negative, dark = (
        palette
    )

    df = results.copy()

    df["plot_p"] = (
        df["phylo_p"]
        .fillna(1.0)
        .clip(lower=1e-300)
    )

    df["minus_log10_p"] = (
        -np.log10(
            df["plot_p"]
        )
    )

    fig, ax = plt.subplots(
        figsize=(8.2, 6.4)
    )

    # Background: all tested KOs
    ax.scatter(
        df["log2_OR"],
        df["minus_log10_p"],
        s=28,
        alpha=0.45,
        c=light,
        edgecolors=dark,
        linewidths=0.35
    )

    # Significant enrichment
    enriched = df[
        df["phylo_significant"]
        & (
            df[
                "prevalence_difference"
            ] > 0
        )
    ]

    ax.scatter(
        enriched["log2_OR"],
        enriched["minus_log10_p"],
        s=48,
        c=positive,
        edgecolors=dark,
        linewidths=0.5,
        label="Enriched with PAPI29"
    )

    # Significant depletion
    depleted = df[
        df["phylo_significant"]
        & (
            df[
                "prevalence_difference"
            ] < 0
        )
    ]

    ax.scatter(
        depleted["log2_OR"],
        depleted["minus_log10_p"],
        s=48,
        c=negative,
        edgecolors=dark,
        linewidths=0.5,
        label="Depleted with PAPI29"
    )

    ax.axvline(
        0,
        color=dark,
        linewidth=0.8,
        alpha=0.5
    )

    # Label strongest adjusted associations
    label_df = (
        df[
            df["phylo_p"].notna()
        ]
        .sort_values(
            "phylo_p"
        )
        .head(12)
    )

    for _, r in label_df.iterrows():

        ax.annotate(
            r["KO"],
            (
                r["log2_OR"],
                r["minus_log10_p"]
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
            color=dark
        )

    ax.set_xlabel(
        "Association effect (log$_2$ odds ratio)"
    )

    ax.set_ylabel(
        "Phylogeny-adjusted −log$_{10}$(P)"
    )

    ax.set_title(
        "LPS/O-antigen associations with PAPI29"
    )

    if (
        len(enriched)
        or len(depleted)
    ):
        ax.legend(
            frameon=False,
            fontsize=8
        )

    ax.spines[
        ["top", "right"]
    ].set_visible(False)

    fig.tight_layout()

    save_figure(
        fig,
        outdir,
        "01_association_scatter"
    )


# ============================================================
# SELECT KOs FOR DISPLAY
# ============================================================

def select_display_kos(
    results,
    top_n
):

    # Prefer phylogenetically significant KOs.
    sig = results[
        results["phylo_significant"]
    ].copy()

    if len(sig):

        sig = sig.sort_values(
            [
                "phylo_FDR",
                "phylo_p"
            ]
        )

        return (
            sig.head(top_n)["KO"]
            .tolist()
        )

    # If none survive FDR, display strongest
    # adjusted signals, clearly without calling
    # them significant.
    usable = results[
        results["phylo_p"].notna()
    ].copy()

    usable = usable.sort_values(
        "phylo_p"
    )

    return (
        usable.head(top_n)["KO"]
        .tolist()
    )


# ============================================================
# FIGURE 2:
# PREVALENCE HEATMAP
# ============================================================

def prevalence_heatmap(
    results,
    kos,
    outdir,
    palette
):

    light, positive, negative, dark = (
        palette
    )

    if not kos:
        return

    df = (
        results
        .set_index("KO")
        .loc[kos]
        .copy()
    )

    values = np.column_stack([
        df[
            "prevalence_PAPI29_neg"
        ].to_numpy(),
        df[
            "prevalence_PAPI29_pos"
        ].to_numpy()
    ])

    # Sort by difference for easier reading
    order = np.argsort(
        df[
            "prevalence_difference"
        ].to_numpy()
    )

    values = values[order]

    labels = (
        df.index.to_numpy()[order]
    )

    fig_height = max(
        4.5,
        0.35 * len(labels) + 1.5
    )

    fig, ax = plt.subplots(
        figsize=(5.8, fig_height)
    )

    # Continuous neutral-to-positive map
    from matplotlib.colors import (
        LinearSegmentedColormap
    )

    cmap = LinearSegmentedColormap.from_list(
        "prevalence",
        [light, positive]
    )

    im = ax.imshow(
        values,
        aspect="auto",
        vmin=0,
        vmax=1,
        cmap=cmap
    )

    ax.set_xticks(
        [0, 1]
    )

    ax.set_xticklabels([
        "PAPI29−",
        "PAPI29+"
    ])

    ax.set_yticks(
        np.arange(
            len(labels)
        )
    )

    ax.set_yticklabels(
        labels,
        fontsize=8
    )

    ax.set_title(
        "Prevalence of selected LPS/O-antigen KOs"
    )

    for i in range(
        values.shape[0]
    ):

        for j in range(
            values.shape[1]
        ):

            value = values[i, j]

            text_colour = (
                "white"
                if value > 0.55
                else dark
            )

            ax.text(
                j,
                i,
                f"{100 * value:.1f}%",
                ha="center",
                va="center",
                fontsize=7,
                color=text_colour
            )

    cbar = fig.colorbar(
        im,
        ax=ax,
        fraction=0.045,
        pad=0.04
    )

    cbar.set_label(
        "Genome prevalence"
    )

    ax.tick_params(
        length=0
    )

    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.tight_layout()

    save_figure(
        fig,
        outdir,
        "02_prevalence_heatmap"
    )


# ============================================================
# TREE ORDER
# ============================================================

def tree_tip_order(
    tree,
    allowed
):

    allowed = set(
        allowed
    )

    order = []

    for terminal in tree.get_terminals():

        acc = accession(
            terminal.name
        )

        if acc in allowed:
            order.append(acc)

    return order


# ============================================================
# FIGURE 3:
# TREE + PAPI29 + KO HEATMAP
# ============================================================

def tree_heatmap(
    tree,
    matrix,
    kos,
    outdir,
    palette
):

    light, positive, negative, dark = (
        palette
    )

    if not kos:
        return

    order = tree_tip_order(
        tree,
        matrix[
            "assembly_accession"
        ]
    )

    if not order:
        return

    indexed = (
        matrix
        .set_index(
            "assembly_accession"
        )
        .loc[order]
    )

    # --------------------------------------------------------
    # Build tree coordinates ourselves so rows align exactly.
    # --------------------------------------------------------

    terminals = [
        t for t in tree.get_terminals()
        if accession(t.name) in set(order)
    ]

    y_positions = {
        acc: i
        for i, acc in enumerate(order)
    }

    # X = distance from root
    depths = tree.depths()

    if not max(
        depths.values()
    ):
        depths = tree.depths(
            unit_branch_lengths=True
        )

    # Determine Y for every clade recursively
    clade_y = {}

    def get_y(clade):

        if clade in clade_y:
            return clade_y[clade]

        if clade.is_terminal():

            acc = accession(
                clade.name
            )

            if acc in y_positions:
                y = y_positions[acc]
            else:
                y = np.nan

        else:

            children = [
                get_y(c)
                for c in clade.clades
            ]

            children = [
                y for y in children
                if np.isfinite(y)
            ]

            if children:
                y = np.mean(children)
            else:
                y = np.nan

        clade_y[clade] = y

        return y

    get_y(
        tree.root
    )

    # --------------------------------------------------------
    # Layout
    # --------------------------------------------------------

    n = len(order)

    fig_height = max(
        8,
        min(
            40,
            n * 0.035
        )
    )

    fig = plt.figure(
        figsize=(
            10 + 0.35 * len(kos),
            fig_height
        )
    )

    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[
            4,
            max(
                2,
                len(kos) * 0.32
            )
        ],
        wspace=0.02
    )

    ax_tree = fig.add_subplot(
        gs[0, 0]
    )

    ax_heat = fig.add_subplot(
        gs[0, 1]
    )

    # --------------------------------------------------------
    # Draw tree
    # --------------------------------------------------------

    def draw_clade(clade):

        x_parent = depths[
            clade
        ]

        valid_children = []

        for child in clade.clades:

            y_child = clade_y.get(
                child,
                np.nan
            )

            if not np.isfinite(
                y_child
            ):
                continue

            valid_children.append(
                child
            )

            x_child = depths[
                child
            ]

            ax_tree.plot(
                [
                    x_parent,
                    x_child
                ],
                [
                    y_child,
                    y_child
                ],
                color=dark,
                linewidth=0.35
            )

            draw_clade(
                child
            )

        if valid_children:

            ys = [
                clade_y[c]
                for c in valid_children
            ]

            ax_tree.plot(
                [
                    x_parent,
                    x_parent
                ],
                [
                    min(ys),
                    max(ys)
                ],
                color=dark,
                linewidth=0.35
            )

    draw_clade(
        tree.root
    )

    ax_tree.set_ylim(
        -1,
        n
    )

    ax_tree.invert_yaxis()

    ax_tree.set_xlabel(
        "Phylogenetic distance"
    )

    ax_tree.set_yticks([])

    ax_tree.spines[
        ["top", "right", "left"]
    ].set_visible(False)

    # --------------------------------------------------------
    # Heatmap
    #
    # First column = PAPI29
    # Remaining columns = selected KOs
    # --------------------------------------------------------

    heat = np.column_stack([
        indexed[
            "PAPI29"
        ].to_numpy(
            dtype=int
        )
    ] + [
        indexed[
            ko
        ].to_numpy(
            dtype=int
        )
        for ko in kos
    ])

    # 0 = light
    # 1 = coloured
    cmap = ListedColormap([
        light,
        positive
    ])

    ax_heat.imshow(
        heat,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=1
    )

    labels = (
        ["PAPI29"]
        + kos
    )

    ax_heat.set_xticks(
        np.arange(
            len(labels)
        )
    )

    ax_heat.set_xticklabels(
        labels,
        rotation=90,
        fontsize=7
    )

    ax_heat.set_yticks([])

    ax_heat.tick_params(
        length=0
    )

    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    fig.suptitle(
        "Phylogenetic distribution of PAPI29 and associated LPS/O-antigen KOs",
        fontsize=11
    )

    save_figure(
        fig,
        outdir,
        "03_tree_KO_heatmap"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    p = argparse.ArgumentParser(
        description=(
            "Test associations between PAPI29 "
            "and LPS/O-antigen KOs while "
            "accounting for phylogenetic "
            "population structure."
        )
    )

    p.add_argument(
        "--tree",
        type=Path,
        required=True,
        help=(
            "GToTree/IQ-TREE Newick tree."
        )
    )

    p.add_argument(
        "--kofam",
        type=Path,
        required=True,
        help=(
            "KOfam binary KO presence matrix "
            "(ko_presence_unambiguous.tsv)."
        )
    )

    p.add_argument(
        "--ipg",
        type=Path,
        required=True,
        help=(
            "PAPI29 IPG distribution table."
        )
    )

    p.add_argument(
        "--cohort",
        type=Path,
        required=True,
        help=(
            "Text file containing the GCF/GCA "
            "accessions defining the cohort."
        )
    )

    p.add_argument(
        "--outdir",
        type=Path,
        required=True,
        help=(
            "Output directory."
        )
    )

    p.add_argument(
        "--palette",
        nargs=4,
        metavar=(
            "LIGHT",
            "POSITIVE",
            "NEGATIVE",
            "DARK"
        ),
        default=[
            "#F5EBDD",
            "#F2765E",
            "#315B8C",
            "#413333"
        ],
        help=(
            "Four plotting colours: "
            "LIGHT POSITIVE NEGATIVE DARK. "
            "Default: "
            "#F5EBDD #F2765E "
            "#315B8C #413333"
        )
    )

    p.add_argument(
        "--phylo-variance",
        type=float,
        default=0.80,
        help=(
            "Fraction of positive phylogenetic "
            "eigenvalue variance retained as "
            "population-structure covariates. "
            "Default: 0.80"
        )
    )

    p.add_argument(
        "--top-n",
        type=int,
        default=20,
        help=(
            "Maximum number of KOs shown in "
            "heatmap/tree figures. "
            "Default: 20"
        )
    )

    args = p.parse_args()

    # --------------------------------------------------------
    # Validate arguments
    # --------------------------------------------------------

    for path in [
        args.tree,
        args.kofam,
        args.ipg,
        args.cohort
    ]:

        if not path.is_file():
            p.error(
                f"Input file does not exist: {path}"
            )

    if not (
        0 < args.phylo_variance <= 1
    ):
        p.error(
            "--phylo-variance must be >0 and <=1."
        )

    if args.top_n < 1:
        p.error(
            "--top-n must be >=1."
        )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Read cohort
    # --------------------------------------------------------

    cohort = read_cohort(
        args.cohort
    )

    cohort_set = set(
        cohort
    )

    print(
        "\n========================================"
    )

    print(
        "PAPI29 × LPS ASSOCIATION ANALYSIS"
    )

    print(
        "========================================"
    )

    print(
        f"Cohort genomes: {len(cohort)}"
    )

    # --------------------------------------------------------
    # Read PAPI29 status
    # --------------------------------------------------------

    papi29_all = read_papi29_ipg(
        args.ipg
    )

    papi29 = (
        papi29_all
        & cohort_set
    )

    print(
        f"PAPI29+ cohort genomes: {len(papi29)}"
    )

    print(
        "PAPI29- cohort genomes:",
        len(cohort) - len(papi29)
    )

    if len(papi29) == 0:
        raise ValueError(
            "No cohort genome was identified "
            "as PAPI29+ from the IPG."
        )

    # --------------------------------------------------------
    # KOfam
    # --------------------------------------------------------

    kofam, ko_columns = (
        read_kofam(
            args.kofam
        )
    )

    kofam_set = set(
        kofam[
            "assembly_accession"
        ]
    )

    missing_kofam = (
        cohort_set
        - kofam_set
    )

    extra_kofam = (
        kofam_set
        - cohort_set
    )

    print(
        f"KOfam genomes: {len(kofam_set)}"
    )

    print(
        f"KOfam KOs tested: {len(ko_columns)}"
    )

    print(
        "Cohort genomes missing from KOfam:",
        len(missing_kofam)
    )

    print(
        "KOfam genomes outside cohort:",
        len(extra_kofam)
    )

    if missing_kofam:

        print(
            "\nFirst missing KOfam genomes:"
        )

        for x in sorted(
            missing_kofam
        )[:20]:
            print(" ", x)

        raise ValueError(
            "Cohort/KOfam mismatch. "
            "Stopping rather than silently "
            "dropping genomes."
        )

    # --------------------------------------------------------
    # Tree
    # --------------------------------------------------------

    print(
        "\nReading tree..."
    )

    tree = Phylo.read(
        args.tree,
        "newick"
    )

    terminals = get_tree_terminals(
        tree
    )

    tree_set = set(
        terminals
    )

    missing_tree = (
        cohort_set
        - tree_set
    )

    extra_tree = (
        tree_set
        - cohort_set
    )

    print(
        f"Tree genomes: {len(tree_set)}"
    )

    print(
        "Cohort genomes missing from tree:",
        len(missing_tree)
    )

    print(
        "Tree genomes outside cohort:",
        len(extra_tree)
    )

    if missing_tree:

        print(
            "\nFirst missing tree genomes:"
        )

        for x in sorted(
            missing_tree
        )[:20]:
            print(" ", x)

        raise ValueError(
            "Cohort/tree mismatch. "
            "Stopping rather than silently "
            "dropping genomes."
        )

    # --------------------------------------------------------
    # Construct exact ordered matrix
    # --------------------------------------------------------

    kofam = (
        kofam
        .set_index(
            "assembly_accession"
        )
        .loc[cohort]
        .reset_index()
    )

    kofam[
        "PAPI29"
    ] = (
        kofam[
            "assembly_accession"
        ]
        .isin(papi29)
        .astype(int)
    )

    # Save exact analysed cohort
    cohort_output = kofam[
        [
            "assembly_accession",
            "PAPI29"
        ]
    ].copy()

    cohort_output.to_csv(
        args.outdir / "cohort.tsv",
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Tree distances
    # --------------------------------------------------------

    D = tree_distance_matrix(
        tree,
        cohort,
        terminals
    )

    np.save(
        args.outdir
        / "patristic_distance_matrix.npy",
        D
    )

    # --------------------------------------------------------
    # Phylogenetic coordinates
    # --------------------------------------------------------

    (
        phylo,
        eigenvalues,
        explained,
        n_axes
    ) = phylogenetic_coordinates(
        D,
        args.phylo_variance
    )

    phylo_df = pd.DataFrame({
        "assembly_accession":
            cohort
    })

    for i in range(
        phylo.shape[1]
    ):

        phylo_df[
            f"PhyloPC{i + 1}"
        ] = phylo[:, i]

    phylo_df.to_csv(
        args.outdir
        / "phylogenetic_coordinates.tsv",
        sep="\t",
        index=False
    )

    eigen_df = pd.DataFrame({
        "axis":
            np.arange(
                1,
                len(eigenvalues) + 1
            ),

        "eigenvalue":
            eigenvalues,

        "variance_fraction":
            explained,

        "cumulative_variance":
            np.cumsum(explained),

        "retained":
            np.arange(
                1,
                len(eigenvalues) + 1
            ) <= n_axes
    })

    eigen_df.to_csv(
        args.outdir
        / "phylogenetic_eigenvalues.tsv",
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Associations
    # --------------------------------------------------------

    results = run_associations(
        kofam,
        ko_columns,
        phylo
    )

    # Sort primarily by adjusted phylogenetic
    # evidence, then naive evidence.
    results = results.sort_values(
        [
            "phylo_FDR",
            "phylo_p",
            "fisher_FDR",
            "fisher_p"
        ],
        na_position="last"
    )

    results.to_csv(
        args.outdir
        / "lps_association_results.tsv",
        sep="\t",
        index=False
    )

    significant = results[
        results[
            "phylo_significant"
        ]
    ].copy()

    significant.to_csv(
        args.outdir
        / "significant_phylogeny_adjusted_associations.tsv",
        sep="\t",
        index=False
    )

    # Also explicitly save things that look
    # significant naively but disappear after
    # population-structure adjustment.
    lineage_dependent = results[
        results["naive_significant"]
        & ~results["phylo_significant"]
    ].copy()

    lineage_dependent.to_csv(
        args.outdir
        / "naive_only_associations.tsv",
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Display KOs
    # --------------------------------------------------------

    display_kos = select_display_kos(
        results,
        args.top_n
    )

    pd.DataFrame({
        "KO": display_kos
    }).to_csv(
        args.outdir
        / "KOs_displayed_in_figures.tsv",
        sep="\t",
        index=False
    )

    # --------------------------------------------------------
    # Figures
    # --------------------------------------------------------

    print(
        "\nDrawing figures...",
        flush=True
    )

    association_scatter(
        results,
        args.outdir,
        args.palette
    )

    prevalence_heatmap(
        results,
        display_kos,
        args.outdir,
        args.palette
    )

    tree_heatmap(
        tree,
        kofam,
        display_kos,
        args.outdir,
        args.palette
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary = {
        "cohort_genomes":
            len(cohort),

        "PAPI29_positive":
            int(
                kofam["PAPI29"].sum()
            ),

        "PAPI29_negative":
            int(
                len(kofam)
                - kofam["PAPI29"].sum()
            ),

        "KOs_tested":
            len(ko_columns),

        "phylogenetic_axes_retained":
            n_axes,

        "phylogenetic_variance_target":
            args.phylo_variance,

        "phylogenetic_variance_represented":
            float(
                np.cumsum(
                    explained
                )[n_axes - 1]
            ),

        "naive_FDR_significant":
            int(
                results[
                    "naive_significant"
                ].sum()
            ),

        "phylogeny_adjusted_FDR_significant":
            int(
                results[
                    "phylo_significant"
                ].sum()
            ),

        "naive_only_significant":
            int(
                (
                    results[
                        "naive_significant"
                    ]
                    & ~results[
                        "phylo_significant"
                    ]
                ).sum()
            )
    }

    pd.Series(
        summary
    ).to_csv(
        args.outdir
        / "analysis_summary.tsv",
        sep="\t",
        header=False
    )

    print(
        "\n========================================"
    )

    print(
        "COMPLETE"
    )

    print(
        "========================================"
    )

    for key, value in summary.items():
        print(
            f"{key}: {value}"
        )

    print(
        "\nResults:"
    )

    print(
        args.outdir
        / "lps_association_results.tsv"
    )

    print(
        "\nPhylogeny-adjusted significant:"
    )

    print(
        args.outdir
        / "significant_phylogeny_adjusted_associations.tsv"
    )

    print(
        "\nFigures:"
    )

    print(
        "  01_association_scatter"
    )

    print(
        "  02_prevalence_heatmap"
    )

    print(
        "  03_tree_KO_heatmap"
    )


if __name__ == "__main__":
    main()
