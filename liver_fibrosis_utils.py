import squidpy as sq
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import os
from crisp import plot as crisp_plot


# Canonical display names for the three interaction-score granularities used
# throughout this project. Any raw column/key name that represents a score
# computed by crisp should be routed through _score_label() so every plot
# uses consistent, human-readable text instead of the internal variable name.
_SCORE_LABELS = {
    "proj_pearson_cell": "Cell-interaction score",
    "proj_spearman_cell": "Cell-interaction score",
    "cell_score": "Cell-interaction score",
    "proj_mean_cell_pearson": "Sample-interaction score",
    "proj_spearman_global": "Sample-interaction score",
    "proj_sim_global": "Sample-interaction score",
    "ratio": "Sample-interaction score",
    "proj_pearson_gene_contrib": "Gene-interaction score",
    "proj_pearson_gene_contrib_per_expr": "Gene-interaction score (per expression)",
    "cs_gene_score": "Gene-interaction score",
}


def _score_label(key):
    """Map a raw score column/key name to its canonical display label."""
    return _SCORE_LABELS.get(key, key)


def apply_plot_style():
    sns.set_theme(style="whitegrid")
    plt.rcParams.update({
        "figure.figsize": (8,4),
        "font.size": 18,
        "axes.titlesize": 24,
        "axes.titleweight": "bold",
        "axes.labelsize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 16,
        "legend.title_fontsize": 17,
        "figure.titlesize": 26,
        "figure.titleweight": "bold",
    })
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"
    plt.rcParams["savefig.edgecolor"] = "white"
    sc.settings.set_figure_params(facecolor="white", fontsize=22, dpi=100)
    return sc.settings.figdir




def plot_cell_type_distribution_by_score_bins(adata_hep, cell_score_types, n_bins, remove_outliers, bin_appear=False):
# Loop over unique samples

    apply_plot_style()
    adata_filt = adata_hep.copy()
    for cell_score in cell_score_types:
        cell_types_of_interest = ["Fibrotic Hep 1", "Fibrotic Hep 2"]
        if remove_outliers:
            score = adata_hep.obs[cell_score]
            low_thr = score.quantile(0.025)
            high_thr = score.quantile(0.975)

            mask = (score >= low_thr) & (score <= high_thr)
            adata_filt = adata_hep[mask].copy()
        print(adata_hep.obs)


        cell_score_bin2 = pd.qcut(
            adata_filt.obs[cell_score],
            q=n_bins,
            duplicates="drop"
        )

        P_ct_given_s_bin2 = pd.crosstab(
            cell_score_bin2,
            adata_hep.obs["Cell_Type_final"],
            normalize="index"
        )

        P_ct_given_s_bin2 = P_ct_given_s_bin2.reindex(
            columns=adata_hep.obs["Cell_Type_final"].unique(),
            fill_value=0
        )

        # ---- Stacked bar ----
        fig, ax = plt.subplots(figsize=(8, 4))

        P_ct_given_s_bin2.plot(kind="bar", stacked=True, ax=ax)

        ax.legend(
            title="Cell Type",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=16,
            title_fontsize=17,
            frameon=False
        )

        if not bin_appear:
            ax.set_xticklabels(
                [f"{i + 1}" for i in range(len(P_ct_given_s_bin2.index))],
                rotation=45
            )
        else:
            ax.set_xticklabels(
                P_ct_given_s_bin2.index.astype(str),
                rotation=45
            )

        ax.set_xlabel("Cell-interaction score bin")
        ax.set_ylabel("P(Cell Type)")
        # ax.set_title(f"P(Cell Type | {_score_label(cell_score)})")

        plt.tight_layout()
        crisp_plot.save_plot(f"{sc.settings.figdir}/cell_type_dist_bin2_{cell_score}.png")
        plt.close()

        # ---- Fibrotic curve ----
        P_selected = (
            P_ct_given_s_bin2.reindex(columns=cell_types_of_interest, fill_value=0)
            .sum(axis=1)
        )

        for bin_name, prob in P_selected.items():
            print(f"Bin: {bin_name} → P({' or '.join(cell_types_of_interest)}) = {prob:.4f}")

        fig, ax = plt.subplots(figsize=(8, 4))

        if not bin_appear:
            x_labels = [f"{i + 1}" for i in range(len(P_selected))]
        else:
            x_labels = P_selected.index.astype(str)

        ax.plot(x_labels, P_selected.values, marker="o")

        ax.set_xlabel("Cell-interaction score bin")
        ax.set_ylabel("P(Fibrotic Hep)")
        # ax.set_title(f"P(Fibrotic Hep) | cell-interaction score bin)")
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, linestyle="--", alpha=0.7)

        plt.tight_layout()
        crisp_plot.save_plot(f"{sc.settings.figdir}/fibrotic_curve_bin2_{cell_score}_ro={remove_outliers}.png")
        plt.close()

        print(cell_score_bin2.value_counts().sort_index())

        # =====================================================
        # --------------------- UMAP --------------------------
        # =====================================================

        sc.pp.neighbors(adata_hep, use_rep="X_pca", n_neighbors=15)
        sc.tl.umap(adata_hep)
        score = adata_hep.obs[cell_score]
        low_thr = score.quantile(0.025)
        high_thr = score.quantile(0.975)

        mask = (score >= low_thr) & (score <= high_thr)
        adata_filt = adata_hep[mask].copy()
        crisp_plot.plot_umap(
            adata_hep,
            save=f"_{cell_score}.png",
            color=[cell_score, "Cell_Type_final"],
            title=[_score_label(cell_score), "Cell Type"],
            vmin = low_thr,
            vmax = high_thr

        )

def plot_umap(labels):
    outdir = apply_plot_style()
    for sample_folder in os.listdir(outdir):

        subdir = os.path.join(outdir, sample_folder)
        if not os.path.isdir(subdir):
            continue
        adata_files = [f for f in os.listdir(subdir) if f.endswith(".h5ad")]

        if not adata_files:
            continue
        import scanpy as sc
        for f in adata_files:
            print(f)
            adata_path = os.path.join(subdir, f)
            adata_tmp = sc.read_h5ad(adata_path)
            sid = sample_folder.replace("sample_", "")
            ratio = adata_tmp.uns.get("ratio", float("nan"))

            print(f"Plotting sample: {sid}, ratio: {ratio:.4f}")

            sc.pp.neighbors(adata_tmp, use_rep="X_pca", n_neighbors=15)
            sc.tl.umap(adata_tmp)

            # ------------------------------------------------------------------
            # label can now be a list
            # example: labels = ["zonation_position", "ALB", "CYP2E1"]
            # ------------------------------------------------------------------

            labels_to_plot = []

            for lab in labels:

                # case 1 → gene
                if lab in adata_tmp.var_names:
                    values = adata_tmp[:, lab].X
                    if hasattr(values, "toarray"):
                        values = values.toarray().flatten()

                    new_name = f"{lab}_expr"
                    adata_tmp.obs[new_name] = values
                    labels_to_plot.append(new_name)

                # case 2 → score / metadata
                else:
                    labels_to_plot.append(lab)

            # ------------------------------------------------------------------
            # compute thresholds for each label separately
            # ------------------------------------------------------------------
            vmin = []
            vmax = []

            for lab in labels_to_plot:
                score = adata_tmp.obs[lab]
                vmin.append(score.quantile(0.025))
                vmax.append(score.quantile(0.975))

            # ------------------------------------------------------------------
            # plot all labels at once
            # ------------------------------------------------------------------
            crisp_plot.plot_umap(
                adata_tmp,
                save="_multi_labels.png",
                color=labels_to_plot + ["Cell_Type_final"],
                title=labels + ["Cell Type"],
                vmin=vmin + [None],
                vmax=vmax + [None]
            )





def plot_label_vs_spatial(label, show_score=False):
    outdir = apply_plot_style()
    # ---------------- Diseased ----------------
    for sample_folder in os.listdir(outdir):

        subdir = os.path.join(outdir, sample_folder)
        if not os.path.isdir(subdir):
            continue
        adata_files = [f for f in os.listdir(subdir) if f.endswith(".h5ad")]

        if not adata_files:
            continue

        adata_path = os.path.join(subdir, adata_files[0])
        adata_tmp = sc.read_h5ad(adata_path)

        sid = sample_folder.replace("sample_", "")
        ratio = adata_tmp.uns.get("proj_spearman_global", float("nan"))

        print(f"Plotting sample: {sid}, ratio: {ratio:.4f}")
        title = f"{sid} (sample-interaction score={ratio:.3f})" if show_score else f"{sid} "
        # add_fibrotic_label(adata_tmp)
        # adata_tmp.write_h5ad(adata_path)


        sq.pl.spatial_scatter(
            adata_tmp,
            color=label,
            colorbar=True,
            shape=None,
            title=title
        )


        fig = plt.gcf()
        ax = fig.axes[0]
        cbar_ax = fig.axes[1]

        # Align colorbar
        cbar_ax.set_position([
            ax.get_position().x1 + 0.02,
            ax.get_position().y0,
            0.03,
            ax.get_position().height
        ])
        if label == "zonation_position":
            label_for_plot = "zonation"
        cbar_ax.set_ylabel(label_for_plot, fontsize=22)

        fig.savefig(
            os.path.join(subdir, f"_{label}_vs_spatial.png"),
            dpi=300,
            bbox_inches='tight'
        )

        plt.close(fig)

def plot_top_genes_spatial_scatter(sample_id, gene_score_key="proj_pearson_gene_contrib", n_top=10):
    """
    Spatial scatter of a single sample, one plot per top-ranked gene (colored by expression).

    Plots are saved under sample_<sample_id>/top_genes_spatial/.
    """
    outdir = apply_plot_style()
    sample_dir = os.path.join(outdir, f"sample_{sample_id}")
    adata_path = os.path.join(sample_dir, "adata.h5ad")

    if not os.path.exists(adata_path):
        print(f"Warning: {adata_path} not found. Skipping.")
        return

    adata_tmp = sc.read_h5ad(adata_path)

    if gene_score_key not in adata_tmp.var.columns:
        print(f"Sample {sample_id}: '{gene_score_key}' not found in adata.var, skipping.")
        return

    subdir = os.path.join(sample_dir, "top_genes_spatial")
    os.makedirs(subdir, exist_ok=True)

    scores = adata_tmp.var[gene_score_key].dropna()
    top_genes = scores.nlargest(n_top)

    for rank, (gene, val) in enumerate(top_genes.items(), start=1):
        print(f"Plotting sample: {sample_id}, rank {rank}: {gene} ({gene_score_key}={val:.4f})")

        sq.pl.spatial_scatter(
            adata_tmp,
            color=gene,
            colorbar=True,
            shape=None,
            # title=f"{sample_id} — {gene}"
            title=f"{gene}",
            cmap="YlGn"

        )

        fig = plt.gcf()
        ax = fig.axes[0]
        cbar_ax = fig.axes[1]

        # Align colorbar
        cbar_ax.set_position([
            ax.get_position().x1 + 0.02,
            ax.get_position().y0,
            0.03,
            ax.get_position().height
        ])

        cbar_ax.set_ylabel("Expression", fontsize=22)

        fig.savefig(
            os.path.join(subdir, f"top{rank}_{gene}_spatial.png"),
            dpi=300,
            bbox_inches='tight'
        )

        plt.close(fig)


def load_global_score(score_type, sample_ids, outdir):
    import numpy as np
    scores = []
    print(sample_ids)
    for sid in sample_ids:
        path = os.path.join(outdir, f"sample_{sid}", f"adata.h5ad")
        if os.path.exists(path):
            adata_tmp = sc.read_h5ad(path)
            print(adata_tmp.n_obs)
            score = adata_tmp.uns[score_type]
            print(score)
            if score is not None:
                scores.append(score)
    return np.array(scores)

def plot_cell_score_box_plots(cell_score):
    outdir = apply_plot_style()
    for sample_folder in os.listdir(outdir):
        subdir = os.path.join(outdir, sample_folder)
        if not os.path.isdir(subdir):
            continue
        adata_files = [f for f in os.listdir(subdir) if f.endswith(".h5ad")]

        if not adata_files:
            continue

        adata_path = os.path.join(subdir, adata_files[0])
        adata_tmp = sc.read_h5ad(adata_path)

        sid = sample_folder.replace("sample_", "")
        ratio = adata_tmp.uns.get("ratio", float("nan"))
        add_fibrotic_label(adata_tmp)
        print(f"Plotting sample: {sid}, ratio: {ratio:.4f}")
        healthy = adata_tmp.obs.loc[~adata_tmp.obs["is_fibrotic"], cell_score]
        diseased = adata_tmp.obs.loc[adata_tmp.obs["is_fibrotic"], cell_score]

        # --- Violin (with inner quartile box). cut=0 clips the KDE to the
        # observed data range so the density does not bleed past the real min/max,
        # which would misrepresent the low tail we specifically want to show. ---
        df_plot = pd.DataFrame({
            "score": pd.concat([healthy, diseased]).to_numpy(dtype=float),
            "group": ["Healthy"] * len(healthy) + ["Fibrotic"] * len(diseased),
        })

        fig, ax = plt.subplots(figsize=(5, 6))
        sns.violinplot(
            data=df_plot, x="group", y="score", order=["Healthy", "Fibrotic"],
            cut=0, inner="quartile", density_norm="width", ax=ax,
        )
        ax.set_ylabel(_score_label(cell_score))
        ax.set_xlabel("")
        # ax.set_title("Cell-interaction score distributions")

        plt.tight_layout()

        plt.savefig(
            os.path.join(subdir, f"_violin_{cell_score}_healthy_vs_fibrotic.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()

        # --- Box-plot version. showfliers=True (unlike the original) so the
        # low-scoring healthy outliers that fall into the fibrotic range are
        # actually drawn instead of being silently suppressed. ---
        fig, ax = plt.subplots(figsize=(5, 6))
        ax.boxplot(
            [healthy, diseased],
            widths=0.6,
            showfliers=True,
            flierprops=dict(marker="d", markersize=3, alpha=0.3,
                            markerfacecolor="black", markeredgecolor="none"),
        )
        ax.set_xticklabels(["Healthy", "Fibrotic"])
        ax.set_ylabel(_score_label(cell_score))
        # ax.set_title("Cell-interaction score distributions")

        plt.tight_layout()
        plt.savefig(
            os.path.join(subdir, f"_boxplot_{cell_score}_healthy_vs_fibrotic.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()

        # --- Separation metrics for the manuscript (printed per sample) ---
        # AUROC = P(random healthy cell scores higher than random fibrotic cell);
        # it is a monotone rescaling of the Mann-Whitney U statistic, and
        # Cliff's delta = 2*AUROC - 1. MW gives the test; AUROC/Cliff's the
        # effect size (p-values are ~0 at these cell counts and not informative).
        import numpy as np
        from scipy.stats import mannwhitneyu
        from sklearn.metrics import roc_auc_score

        s_nf = healthy.to_numpy(dtype=float)
        s_fib = diseased.to_numpy(dtype=float)
        if len(s_nf) and len(s_fib):
            _, mw_p = mannwhitneyu(s_nf, s_fib, alternative="greater")
            auroc = roc_auc_score(
                [1] * len(s_nf) + [0] * len(s_fib),
                np.concatenate([s_nf, s_fib]),
            )
            cliffs = 2 * auroc - 1
            tau = float(np.median(s_fib))               # fibrotic median
            frac_nf_low = float((s_nf < tau).mean())     # healthy overlap into fibrotic range
            nf_med, nf_q1, nf_q3 = np.percentile(s_nf, [50, 25, 75])
            fib_med, fib_q1, fib_q3 = np.percentile(s_fib, [50, 25, 75])
            print(f"  [metrics] AUROC={auroc:.3f}  Cliff_delta={cliffs:.3f}  "
                  f"MW_p={mw_p:.1e}  tau(fib median)={tau:.3f}  "
                  f"healthy<tau={frac_nf_low*100:.1f}%")
            print(f"  [medians/IQR] healthy: median={nf_med:.3f} "
                  f"IQR=[{nf_q1:.3f}, {nf_q3:.3f}]  |  "
                  f"fibrotic: median={fib_med:.3f} IQR=[{fib_q1:.3f}, {fib_q3:.3f}]")

        print("num of fibrotic cells:", adata_tmp.obs["is_fibrotic"].sum())
        print("num of healthy cells:", healthy.shape[0])


def plot_global_scores(healthy_samples, diseased_samples, global_score_name="ratio"):
    import numpy as np
    import matplotlib.pyplot as plt
    outdir = apply_plot_style()
    print(outdir)

    # -------- Healthy --------
    healthy = load_global_score(global_score_name,healthy_samples, outdir)
    unhealthy = load_global_score(global_score_name, diseased_samples, outdir)
    print(healthy, unhealthy)
    fig, axes = plt.subplots(1, 2, figsize=(8, 5), sharey=True)
    plt.subplots_adjust(wspace=0.5)
    # --- Healthy ---
    axes[0].boxplot(
        healthy,
        widths=0.5,
        showfliers=False
    )

    axes[0].scatter(
        np.ones_like(healthy),
        healthy,
        zorder=3
    )


    # --- Unhealthy ---
    axes[1].boxplot(
        unhealthy,
        widths=0.5,
        showfliers=False
    )

    axes[1].scatter(
        np.ones_like(unhealthy),
        unhealthy,
        zorder=3
    )



    axes[0].set_ylabel("Sample-interaction score")
    axes[0].set_title("Healthy samples")
    axes[1].set_title("Diseased samples")

    axes[0].set_xticks([])
    axes[1].set_xticks([])

    plt.savefig(
        f"{outdir}/healthy_vs_unhealthy_livers.png",
        dpi=300,
        bbox_inches="tight"
    )

def plot_score_vs_fibrotic_fraction(healthy_samples, diseased_samples, global_score_name="proj_spearman_global"):
    import numpy as np
    import matplotlib.pyplot as plt
    outdir = apply_plot_style()

    records = []
    for sid in healthy_samples:
        path = os.path.join(outdir, f"sample_{sid}", "adata.h5ad")
        if not os.path.exists(path):
            continue
        adata_tmp = sc.read_h5ad(path)
        add_fibrotic_label(adata_tmp)
        score = adata_tmp.uns.get(global_score_name)
        if score is None:
            continue
        fibrotic_pct = adata_tmp.obs["is_fibrotic"].mean() * 100
        records.append({"sample": sid, "score": float(score), "fibrotic_pct": fibrotic_pct, "condition": "healthy"})

    for sid in diseased_samples:
        path = os.path.join(outdir, f"sample_{sid}", "adata.h5ad")
        if not os.path.exists(path):
            continue
        adata_tmp = sc.read_h5ad(path)
        add_fibrotic_label(adata_tmp)
        score = adata_tmp.uns.get(global_score_name)
        if score is None:
            continue
        fibrotic_pct = adata_tmp.obs["is_fibrotic"].mean() * 100
        records.append({"sample": sid, "score": float(score), "fibrotic_pct": fibrotic_pct, "condition": "diseased"})

    df = pd.DataFrame(records)

    colors = {"healthy": "#2ca02c", "diseased": "#d62728"}
    fig, ax = plt.subplots(figsize=(7, 5))

    for condition, group in df.groupby("condition"):
        ax.scatter(
            group["fibrotic_pct"],
            group["score"],
            color=colors[condition],
            label=condition.capitalize(),
            s=80,
            zorder=3,
        )
        for _, row in group.iterrows():
            ax.annotate(
                row["sample"],
                (row["fibrotic_pct"], row["score"]),
                fontsize=12,
                ha="left",
                va="bottom",
                xytext=(4, 4),
                textcoords="offset points",
            )

    ax.set_xlabel("% Fibrotic Cells")
    ax.set_ylabel(_score_label(global_score_name))
    ax.set_title("Sample-interaction score vs fibrotic cell fraction")
    ax.legend(title="Condition")
    ax.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    crisp_plot.save_plot(f"{outdir}/score_vs_fibrotic_fraction_{global_score_name}.png")
    plt.show()


def add_fibrotic_label(adata):
    fibrotic_hep = ["Fibrotic Hep 1", "Fibrotic Hep 2"]
    adata.obs["is_fibrotic"] = adata.obs["Cell_Type_final"].isin(fibrotic_hep).astype(int)
    adata.obs["is_fibrotic"] = adata.obs["is_fibrotic"].astype(bool)


def _build_spatial_clusters(adata, res=1):
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=15, radius=None)
    sc.tl.leiden(
        adata,
        adjacency=adata.obsp["spatial_connectivities"],
        key_added="spatial_cluster",
        resolution=res,
    )


def _cluster_aggregates(adata, cell_score, gene):
    cluster_means = adata.obs.groupby("spatial_cluster")[cell_score].mean()
    cluster_fibrosis = adata.obs.groupby("spatial_cluster")["is_fibrotic"].mean()
    gene_means = (
        adata.to_df()[[gene]]
        .groupby(adata.obs["spatial_cluster"])
        .mean()
        .squeeze()
    )
    return cluster_means, cluster_fibrosis, gene_means


def _plot_spatial_maps(adata, cluster_means, cluster_fibrosis, cell_score, save_path):
    cluster_scaled = (cluster_means - cluster_means.min()) / (cluster_means.max() - cluster_means.min())
    adata.obs["cluster_scaled"] = adata.obs["spatial_cluster"].map(cluster_scaled).astype(float)
    adata.obs["cluster_fibrotic_prop"] = adata.obs["spatial_cluster"].map(cluster_fibrosis).astype(float)
    sq.pl.spatial_scatter(
        adata,
        color=["cluster_scaled", "cluster_fibrotic_prop"],
        cmap="rocket",
        colorbar=True,
        shape=None,
        ncols=2,
        figsize=(6, 6),
        title=["Mean cell-interaction score per cluster", "Fibrotic proportion per cluster"],
    )
    for ax in plt.gcf().axes:
        ax.title.set_fontsize(24)
        ax.title.set_fontweight("bold")
        ax.tick_params(labelsize=18)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def _scatter_with_spearman(x, y, xlabel, ylabel, title, save_path, ylim=None):
    from scipy.stats import spearmanr
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, alpha=0.7)
    corr, _ = spearmanr(x, y)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} (ρ={corr:.2f})")
    if ylim is not None:
        ax.set_ylim(*ylim)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_mean_cell_score_and_gene_in_spatial_clusters(cell_score, gene, show_ratio=False, remove_outliers=False):
    outdir = apply_plot_style()
    res = 1

    for sample_folder in os.listdir(outdir):
        subdir = os.path.join(outdir, sample_folder)
        if not os.path.isdir(subdir):
            continue
        adata_files = [f for f in os.listdir(subdir) if f.endswith(".h5ad")]
        if not adata_files:
            continue

        adata_tmp = sc.read_h5ad(os.path.join(subdir, adata_files[0]))
        sid = sample_folder.replace("sample_", "")
        ratio = adata_tmp.uns.get("ratio", float("nan"))
        print(f"Plotting sample: {sid}, ratio: {ratio:.4f}")

        add_fibrotic_label(adata_tmp)
        if remove_outliers:
            high_thr = adata_tmp.obs[cell_score].quantile(0.995)
            adata_tmp = adata_tmp[adata_tmp.obs[cell_score] <= high_thr].copy()

        _build_spatial_clusters(adata_tmp, res=res)
        cluster_means, cluster_fibrosis, gene_means = _cluster_aggregates(adata_tmp, cell_score, gene)

        _plot_spatial_maps(
            adata_tmp, cluster_means, cluster_fibrosis, cell_score,
            os.path.join(subdir, f"_spatial_{cell_score}_and_fibrosis_{res}.png"),
        )
        _scatter_with_spearman(
            cluster_means, cluster_fibrosis,
            "Mean cell-interaction score per cluster", "Fibrotic proportion per cluster",
            "Cluster cell-interaction score vs fibrosis",
            os.path.join(subdir, f"_cluster_{cell_score}_vs_fibrosis_{res}.png"),
            ylim=(0, 1),
        )
        _scatter_with_spearman(
            cluster_fibrosis, gene_means,
            "Fibrotic Proportion per Cluster", f"{gene} per Cluster",
            f"Fibrotic Proportion vs {gene}",
            os.path.join(subdir, f"_cluster_fibrosis_vs_{gene}_{res}.png"),
        )



def print_top_genes_by_score(gene_score_key="proj_spearman_gene", n_top=10):
    outdir = apply_plot_style()
    for sample_folder in os.listdir(outdir):
        subdir = os.path.join(outdir, sample_folder)
        if not os.path.isdir(subdir):
            continue
        adata_files = [f for f in os.listdir(subdir) if f.endswith(".h5ad")]

        if not adata_files:
            continue

        adata_path = os.path.join(subdir, adata_files[0])
        adata_tmp = sc.read_h5ad(adata_path)

        sid = sample_folder.replace("sample_", "")
        ratio = adata_tmp.uns.get("ratio", float("nan"))

        if gene_score_key not in adata_tmp.var.columns:
            print(f"Sample {sid}: '{gene_score_key}' not found in adata.var, skipping.")
            continue

        scores = adata_tmp.var[gene_score_key].dropna()
        top_genes = scores.nlargest(n_top)

        print(f"\nSample {sid} (ratio={ratio:.4f}) — top {n_top} genes by {gene_score_key}:")
        for gene, val in top_genes.items():
            print(f"  {gene}: {val:.4f}")



def plot_transition_genes(cell_score=None, n_top=20, early_transition=True,
                          only_healthy_cells=False, gene_score_key=None,
                          corr_type="pearson"):
    """
    For each sample scatter log2FC (fibrotic/non-fibrotic) vs a gene score.

    Gene score priority:
      gene_score_key  → read directly from adata.var[gene_score_key]
      cell_score      → correlate every gene with that cell-level score
                        (Spearman or Pearson, healthy or all cells via only_healthy_cells)
      neither         → raises ValueError

    early_transition=True  → rank by |score| / (|FC| + 0.5)
    early_transition=False → rank by |score| × |FC|
    """
    if gene_score_key is None and cell_score is None:
        raise ValueError("Provide either gene_score_key or cell_score.")

    import numpy as np
    from scipy.stats import rankdata, spearmanr, pearsonr
    outdir = apply_plot_style()

    for sample_folder in os.listdir(outdir):
        subdir = os.path.join(outdir, sample_folder)
        if not os.path.isdir(subdir):
            continue
        adata_files = [f for f in os.listdir(subdir) if f.endswith(".h5ad")]
        if not adata_files:
            continue

        adata_tmp = sc.read_h5ad(os.path.join(subdir, adata_files[0]))
        sid = sample_folder.replace("sample_", "")

        add_fibrotic_label(adata_tmp)
        fibrotic_mask = adata_tmp.obs["is_fibrotic"].values

        if fibrotic_mask.sum() == 0:
            print(f"Sample {sid}: no fibrotic cells, skipping.")
            continue

        # --- gene expression matrices ---
        X_all = adata_tmp.X.toarray() if hasattr(adata_tmp.X, "toarray") else np.array(adata_tmp.X)
        X_nf  = X_all[~fibrotic_mask]
        X_fib = X_all[fibrotic_mask]

        # --- gene y-axis score ---
        if gene_score_key is not None:
            if gene_score_key not in adata_tmp.var.columns:
                print(f"Sample {sid}: '{gene_score_key}' not found in adata.var, skipping.")
                continue
            gene_scores = adata_tmp.var[gene_score_key].values
            y_label = _score_label(gene_score_key)
            score_key_tag = gene_score_key
            fname_tag = f"geneScore={gene_score_key}"
        else:
            if cell_score not in adata_tmp.obs.columns:
                print(f"Sample {sid}: '{cell_score}' not found in adata.obs, skipping.")
                continue
            score = adata_tmp.obs[cell_score].values
            X_corr     = X_nf               if only_healthy_cells else X_all
            score_corr = score[~fibrotic_mask] if only_healthy_cells else score
            cells_tag  = "healthy cells only" if only_healthy_cells else "all cells"

            if corr_type == "pearson":
                score_c = score_corr - score_corr.mean()
                X_c = X_corr - X_corr.mean(axis=0)
                gene_scores = (X_c * score_c[:, None]).sum(axis=0) / (
                    np.linalg.norm(X_c, axis=0) * (np.linalg.norm(score_c) + 1e-12) + 1e-12)
                y_label = f"Gene correlation with cell-interaction score\n({cells_tag})"
            else:
                rank_score = rankdata(score_corr).astype(float)
                rank_score -= rank_score.mean()
                rank_X = np.apply_along_axis(rankdata, 0, X_corr).astype(float)
                rank_X -= rank_X.mean(axis=0)
                gene_scores = (rank_X * rank_score[:, None]).sum(axis=0) / (
                    np.linalg.norm(rank_X, axis=0) * (np.linalg.norm(rank_score) + 1e-12) + 1e-12)
                y_label = f"Gene correlation with cell-interaction score\n({cells_tag})"

            hc_tag = "healthyCells" if only_healthy_cells else "allCells"
            score_key_tag = cell_score
            fname_tag = f"cellScore={cell_score}_{corr_type}_{hc_tag}"

        # --- log2 fold change: fibrotic vs non-fibrotic ---
        mean_fib = X_fib.mean(axis=0) + 1e-6
        mean_nf  = X_nf.mean(axis=0)  + 1e-6
        log2fc   = np.log2(mean_fib / mean_nf)

        results = pd.DataFrame({
            "gene":                  adata_tmp.var_names,
            "gene_score":            gene_scores,
            "log2fc_fibrotic_vs_nf": log2fc,
        }).set_index("gene")

        if early_transition:
            results["ranking_score"] = (
                results["gene_score"].abs() / (results["log2fc_fibrotic_vs_nf"].abs() + 0.5)
            )
            label, suffix = "early-transition", "early"
        else:
            results["ranking_score"] = (
                results["gene_score"].abs() * results["log2fc_fibrotic_vs_nf"].abs()
            )
            label, suffix = "anti-diagonal (DE + continuous)", "antidiag"

        results = results.sort_values("ranking_score", ascending=False)
        results.to_csv(os.path.join(subdir, f"transition_genes_{suffix}_{score_key_tag}.csv"))

        # --- scatter plot ---
        x_vals = results["log2fc_fibrotic_vs_nf"].values
        y_vals = results["gene_score"].values
        corr_val, corr_pval = spearmanr(x_vals, y_vals)
        pear_val, pear_pval = pearsonr(x_vals, y_vals)
        print(f"Sample {sid} [{cells_tag if gene_score_key is None else gene_score_key}]: "
              f"Spearman(FC, gene-score) = {corr_val:.3f} (p={corr_pval:.1e})  |  "
              f"Pearson = {pear_val:.3f} (p={pear_pval:.1e})")
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(x_vals, y_vals, s=20, alpha=0.4, color="steelblue")

        # fit line
        slope, intercept = np.polyfit(x_vals, y_vals, 1)
        x_fit = np.linspace(x_vals.min(), x_vals.max(), 100)
        ax.plot(x_fit, slope * x_fit + intercept, color="firebrick", linewidth=1.5)

        ax.axhline(0, color="black", linewidth=0.7, linestyle="--")
        ax.axvline(0, color="black", linewidth=0.7, linestyle="--")
        ax.set_xlabel("log2 FC  (fibrotic / non-fibrotic)")
        ax.set_ylabel(y_label)
        # ax.set_title(f"Sample {sid} — transition gene landscape\n"
        #              f"Spearman(FC, gene-interaction score) = {corr_val:.3f}")
        plt.tight_layout()
        fig.savefig(os.path.join(subdir, f"transition_genes_{suffix}_{fname_tag}.png"),
                    dpi=300, bbox_inches="tight")
        plt.show()

        print(f"\nSample {sid} — all genes ranked by {label}:")
        print(results[["gene_score", "log2fc_fibrotic_vs_nf", "ranking_score"]].to_string())

