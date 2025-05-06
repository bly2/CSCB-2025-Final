import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
from biomart import BiomartServer
from io import StringIO
import anndata as ad
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from scipy.spatial.distance import cdist
from hmmlearn import hmm
import warnings

# Functions
def fetch_positions(adata):
    # Connect to Ensembl Biomart server
    server = BiomartServer("http://grch37.ensembl.org/biomart")
    dataset = server.datasets['hsapiens_gene_ensembl']

    # Query gene names for only missing gene positions
    no_positions = adata[:, adata.var[['start', 'end']].isna().any(axis=1)]
    with_positions = adata[:, ~adata.var[['start', 'end']].isna().any(axis=1)]
    response = dataset.search({
        'filters':{'ensembl_gene_id':list(no_positions.var['gene_ids'])},
        'attributes':['ensembl_gene_id','chromosome_name','start_position','end_position','strand']
    })

    # Convert response to DataFrame and merge with adata.var if response is successful
    if response.status_code == 200:
        print("Request successful!")
        gene_annotations_df = pd.read_csv(StringIO(response.text),sep='\t',header=None)
        gene_annotations_df.columns = ['gene_ids','chromosome','start','end','strand']
    else:
        print(f"Request failed with status code: {response.status_code}")
        print(response.text)

    # Isolate fetched genes from BioMart in no_positions adata
    fetched_positions = no_positions[:, no_positions.var['gene_ids'].isin(gene_annotations_df['gene_ids'])].copy()

    # Sort fetched genes based on ensembl gene IDs
    fetched_positions = fetched_positions[:, fetched_positions.var['gene_ids'].argsort()].copy()

    # Add the fetched gene positions to the adata
    fetched_positions.var['chromosome'] = gene_annotations_df['chromosome'].values
    fetched_positions.var['start'] = gene_annotations_df['start'].values
    fetched_positions.var['end'] = gene_annotations_df['end'].values
    fetched_positions.var['strand'] = gene_annotations_df['strand'].values

    # Concatenate fetched genes with isolated genes already with positions
    adClean = ad.concat([with_positions, fetched_positions], axis=1)

    # Include obs into the cleaned adata
    adClean.obs = with_positions.obs.copy()

    return adClean

def standardize_chromosomes(adata):
    
    adata1 = adata.copy()

    # Add 'chr' prefix to chromosome names
    adata1.var['chromosome'] = 'chr' + adata1.var['chromosome'].astype(str)

    # Define standard chromosome names with 'chr' prefix
    standard_chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY', 'chrMT']

    # Filter adata to include only genes on standard chromosomes
    adata1 = adata1[:, adata1.var['chromosome'].isin(standard_chromosomes)].copy()

    return adata1

def qc(adata,
    mt_threshold_pct=20,
    min_genes=200,
    max_counts=50000,
    min_cells=3):

    adata1 = adata.copy()

    # Find MT genes
    adata1.var['mt'] = adata1.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata1, qc_vars=['mt'],
                            percent_top=None,
                            log1p=False,
                            inplace=True)

    # Filter out cells based on MT genes
    adClean = adata1[adata1.obs['pct_counts_mt']<mt_threshold_pct,:].copy()

    # Filter out cells based on number of genes expressed
    sc.pp.filter_cells(adClean, min_genes=min_genes)

    # Filter out cells based on total counts
    sc.pp.filter_cells(adClean, max_counts=max_counts)

    # Filter out genes expressed in few cells
    sc.pp.filter_genes(adClean, min_cells=min_cells)

    return adClean

def log_freeman_tukey_transform(expr_mat):
    """
    Freeman-Tukey variance stabilizing transformation:
    log2(sqrt(x) + sqrt(x+1))
    """
    return np.log2(np.sqrt(expr_mat) + np.sqrt(expr_mat + 1))

def identify_diploid_cells_high_precision(adata,
                                           window=10,
                                           n_components=3,
                                           n_pcs=20,
                                           primary_trim_percentile=100,
                                           secondary_trim_percentile=30):
    """
    Identify diploid cells by:
    - Selecting lowest-variance GMM cluster (primary)
    - Adding central portion of second-lowest variance cluster (secondary)
    - Based on PCA + GMM and genomic smoothing
    """
    expr_raw = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    expr_df = pd.DataFrame(expr_raw, index=adata.obs_names, columns=adata.var_names)
    gene_order = adata.var.sort_values(['chromosome', 'start']).index
    expr_df = expr_df[gene_order]

    # Smoothing
    expr_smooth = expr_df.T.rolling(window=window, min_periods=1, center=True).mean().dropna().T
    expr_smooth = expr_smooth.loc[:, ~expr_smooth.columns.duplicated()]

    # PCA
    pca = PCA(n_components=n_pcs)
    pcs = pca.fit_transform(expr_smooth)

    # GMM Clustering
    gmm = GaussianMixture(n_components=n_components, covariance_type='full', random_state=0)
    labels = gmm.fit_predict(pcs)

    # Rank clusters by composite score = variance × size
    cluster_variances = [expr_smooth.iloc[labels == i].var(axis=1).mean() for i in range(n_components)]
    cluster_sizes = [np.sum(labels == i) for i in range(n_components)]
    composite_score = [v * s for v, s in zip(cluster_variances, cluster_sizes)]
    sorted_clusters = np.argsort(composite_score)

    # --- PRIMARY cluster (lowest variance): take central X%
    primary_cluster = sorted_clusters[0]
    primary_idx = np.where(labels == primary_cluster)[0]
    primary_pcs = pcs[primary_idx]
    center1 = primary_pcs.mean(axis=0).reshape(1, -1)
    dist1 = cdist(primary_pcs, center1).flatten()
    thresh1 = np.percentile(dist1, primary_trim_percentile)
    keep_primary = primary_idx[dist1 <= thresh1]

    # --- SECONDARY cluster (2nd lowest variance): take tighter X%
    secondary_cluster = sorted_clusters[1]
    secondary_idx = np.where(labels == secondary_cluster)[0]
    secondary_pcs = pcs[secondary_idx]
    center2 = secondary_pcs.mean(axis=0).reshape(1, -1)
    dist2 = cdist(secondary_pcs, center2).flatten()
    thresh2 = np.percentile(dist2, secondary_trim_percentile)
    keep_secondary = secondary_idx[dist2 <= thresh2]

    # --- Combine both
    confident_indices = np.concatenate([keep_primary, keep_secondary])
    diploid_pred = np.zeros(pcs.shape[0], dtype=bool)
    diploid_pred[confident_indices] = True
    diploid_pred = ['diploid' if x else 'aneuploid' for x in diploid_pred]

    return diploid_pred, labels

def evaluate_predictions(adata, diploid_pred):
    """
    Compare predicted diploids vs. simulated CNV ground truth.
    Returns: precision, recall, F1 score.
    """
    true_diploid = adata.obs['simulated_cnvs'].astype(str).replace(['', 'nan', 'NaN'], np.nan).isna()
    tp = np.sum((diploid_pred == True) & (true_diploid == True))
    fp = np.sum((diploid_pred == True) & (true_diploid == False))
    fn = np.sum((diploid_pred == False) & (true_diploid == True))
    precision = tp / (tp + fp) if (tp + fp) else 0
    recall = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0
    return precision, recall, f1

def downsample(adata,n_cells):
    if adata.n_obs <= n_cells:
        return adata
    else:
        return adata[np.random.choice(adata.obs_names, n_cells, replace=False), :].copy()
    
def cnv_plots(adata,annotation='cell_type'):
    cnv.tl.pca(adata)
    cnv.pp.neighbors(adata)
    cnv.tl.leiden(adata)

    sc.tl.dendrogram(adata, groupby='cnv_leiden')

    cnv.tl.umap(adata)
    cnv.tl.cnv_score(adata)


    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
    ax4.axis("off")
    cnv.pl.umap(
        adata,
        color="cnv_leiden",
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,
    )
    cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(adata, color="cell_type", ax=ax3)

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"wspace": 0.5})
    ax4.axis("off")
    sc.pl.umap(adata, color="cnv_leiden", ax=ax1, show=False)
    sc.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
    sc.pl.umap(adata, color="cell_type", ax=ax3)

def extract_cnv_info(cnv_annotation):
    if cnv_annotation:
        parts = cnv_annotation.split(" ")
        if len(parts) > 1:
            cnv_chr = parts[0].split(":")[0]
            cnv_start = parts[0].split(':')[1].split('-')[0]
            cnv_end = parts[0].split(':')[1].split('-')[1]
            cnv_type = parts[-1].strip("()")
            return cnv_chr,cnv_start,cnv_end,cnv_type
    return None, None