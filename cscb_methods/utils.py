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
import itertools
import warnings

# Functions

def merge_gene_positions(query_adata, ref_adata, gene_colname="gene_ids"):
    '''
    Update gene positions (chromosome, start, end, strand) in the un-annotated dataset with gene positions 
    from an annotated dataset (ref_ad). 
    Also input the name of column containing gene names/ids in gene_colname. It should be the same in both dfs.
    Updates the query_adata, does not return a new object.
    '''
    # get all positions from reference dataset
    ref_pos = ref_adata.var.loc[:, [gene_colname, "chromosome", "start", "end", "strand"]]

    # merge
    pos_overlap = pd.merge(query_adata.var, ref_pos, on=gene_colname, how="left")

    # assign to query datasets .var object
    query_adata.var = pos_overlap

    
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

    # Sort adata by chromosome and start position
    adClean = adClean[:,adClean.var.sort_values(by=['chromosome','start']).index].copy()

    return adClean


def fetch_positions_new(adata, batch_size=200):
    """
    Annotate genes in `adata` with chromosome position info from Ensembl BioMart (GRCh37).
    Fills missing ['chromosome', 'start', 'end', 'strand'] in `.var`, in batches.
    """

    # Connect to Ensembl Biomart server
    server = BiomartServer("http://grch37.ensembl.org/biomart")
    dataset = server.datasets['hsapiens_gene_ensembl']

    # Separate genes with and without positions
    no_positions = adata[:, adata.var[['start', 'end']].isna().any(axis=1)].copy()
    with_positions = adata[:, ~adata.var[['start', 'end']].isna().any(axis=1)].copy()

    # Get gene_ids to query
    gene_ids = no_positions.var['gene_ids'].dropna().unique().tolist()

    # Fetch annotations in batches
    fetched = []
    for i in range(0, len(gene_ids), batch_size):
        batch = gene_ids[i:i + batch_size]
        try:
            response = dataset.search({
                'filters': {'ensembl_gene_id': batch},
                'attributes': ['ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position', 'strand']
            })

            df = pd.read_csv(StringIO(response.text), sep='\t', header=None)
            df.columns = ['gene_ids', 'chromosome', 'start', 'end', 'strand']
            fetched.append(df)
        except Exception as e:
            print(f"Batch {i} failed: {e}")

    # Combine and merge fetched results
    if not fetched:
        print("No annotations fetched.")
        return adata

    gene_annotations_df = pd.concat(fetched, ignore_index=True)

    # Filter and sort no_position genes found in fetched set
    is_fetched = no_positions.var['gene_ids'].isin(gene_annotations_df['gene_ids'])
    fetched_positions = no_positions[:, is_fetched].copy()
    sorted_idx = fetched_positions.var['gene_ids'].argsort()
    fetched_positions = fetched_positions[:, sorted_idx].copy()

    # Map fetched annotations into .var
    gene_annotations_df = gene_annotations_df.set_index(fetched_positions.var.index)
    for col in ['chromosome', 'start', 'end', 'strand']:
        fetched_positions.var[col] = gene_annotations_df[col].values

    # Combine with already-positioned genes
    adClean = ad.concat([with_positions, fetched_positions], axis=1)
    adClean.obs = adata.obs.copy()

    # Optional: sort by chromosome + start
    if 'chromosome' in adClean.var.columns and 'start' in adClean.var.columns:
        adClean = adClean[:, adClean.var.sort_values(['chromosome', 'start']).index].copy()

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
    
def plot_aneuploid_cnv_clusters(adata,diploid_annotation='predicted_diploid'):
    adata_aneuploid = adata[adata.obs['predicted_diploid']=='aneuploid']

    cnv.tl.pca(adata_aneuploid)
    cnv.pp.neighbors(adata_aneuploid)
    cnv.tl.leiden(adata_aneuploid)

    sc.tl.dendrogram(adata_aneuploid, groupby='cnv_leiden')

    cnv.tl.umap(adata_aneuploid)
    cnv.tl.cnv_score(adata_aneuploid)


    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
    ax4.axis("off")
    cnv.pl.umap(
        adata_aneuploid,
        color="cnv_leiden",
        legend_loc="on data",
        legend_fontoutline=2,
        ax=ax1,
        show=False,
    )
    cnv.pl.umap(adata_aneuploid, color="cnv_score", ax=ax2, show=False)
    cnv.pl.umap(adata_aneuploid, color="cell_type", ax=ax3)

def i3_hmm_infercnv(adata,cell_type,cell_annotation='cell_type',diploid_annotation='predicted_diploid',logFC_threshold=0.5,plots=True):
    """ Our main CNV inference approach, utilizing a 3-state Hidden Markov Models like InferCNV to detect the CNV state of a cell (deletion, neutral, amplification) and extract genomic region information and CNV type.

    Args:
        adata (AnnData): AnnData with annotated .obs columns for cell type and diploid/aneuploid cells.
        cell_type (_type_): _description_
        cell_annotation (str, optional): _description_. Defaults to 'cell_type'.
        diploid_annotation (str, optional): _description_. Defaults to 'predicted_diploid'.
        logFC_threshold (float, optional): _description_. Defaults to 0.5.
        plots (bool, optional): _description_. Defaults to True.

    Returns:
        _type_: _description_
    """    
    # Approach similar to 3-state Hidden Markov Model used by the R version of InferCNV to detect genomic regions

    # Transition matrix (3x3) for states: Deletion (0.5), Neutral (1), Amplification (1.5)
    # The transition matrix defines the probability of moving from one state to another. Here, the states are:
    # 0 (Deletion), 1 (Neutral), and 2 (Amplification). The values indicate the transition probabilities between these states.

    transition_matrix = np.array([
        [1 - 5 * 0.000001, 0.000001, 0.000001],  # Transition probabilities for state 0 (Deletion)
        [0.000001, 1 - 5 * 0.000001, 0.000001],  # Transition probabilities for state 1 (Neutral)
        [0.000001, 0.000001, 1 - 5 * 0.000001]   # Transition probabilities for state 2 (Amplification)
    ])

    # Mean expression values for each state (Deletion, Neutral, Amplification)
    means = np.array([0.5, 1.0, 1.5])

    # Covariance values for each state, ensuring they are at least 1e-6 to avoid numerical issues
    covariances = np.maximum(np.array([[0.1], [0.1], [0.1]]), 1e-6)

    # Adata subsetting
    adata_cells = adata[adata.obs[cell_annotation] == cell_type ] 
    diploid_mask = adata_cells.obs['predicted_diploid'] == 'diploid'
    aneuploid_mask = adata_cells.obs['predicted_diploid'] == 'aneuploid'

    exp_all = adata_cells.X.toarray()

    # Compute the diploid reference profile by averaging the expression of diploid cells
    ref_profile = exp_all[diploid_mask].mean(axis=0)

    # Get expression data for aneuploid cells over the background reference
    expX = np.log1p(exp_all[aneuploid_mask] / (ref_profile + 1e-6))

    # Initialize the Hidden Markov Model (HMM) with 3 states (Deletion, Neutral, Amplification)
    model = hmm.GaussianHMM(n_components=3, covariance_type="diag", init_params="")

    # Initialize the starting probabilities for each state. Here, all states have equal probability to start with.
    model.startprob_ = np.array([1/3, 1/3, 1/3])

    # Set the transition matrix and emission probabilities (means and covariances) for the model
    model.transmat_ = transition_matrix
    model.means_ = means.reshape(-1, 1)  # Reshape means to match HMM expectations
    model.covars_ = covariances  # Set covariances for the states

    # Prepare lists to store the predicted CNV levels and states for each cell
    cnv_states_all = []
    predicted_cnv_all = []

    # Define CNV levels based on the hidden states: 0 -> Deletion, 1 -> Neutral, 2 -> Amplification
    cnv_levels = {0: 0.5, 1: 1.0, 2: 1.5}

    # Compute the mean expression across genes (used to filter out low-mean genes)
    gene_means = np.mean(expX, axis=0)

    # Select genes that have a mean expression greater than LogFC threshold
    genes_to_keep = gene_means > 0.5
    expX = expX[:, genes_to_keep]

    # Recompute gene means after filtering, add a small epsilon to avoid division by zero
    gene_means = np.mean(expX, axis=0) + 1e-6

    # Normalize the expression data by dividing by the gene means
    expX = expX / gene_means

    # adata to keep track of gene order
    adata_filtered = adata_cells[aneuploid_mask].copy()
    adata_filtered = adata_filtered[:,genes_to_keep]

    # Iterate through each cell's expression data and predict the CNV states using the HMM model
    for cell_expr in expX:
        cell_expr = cell_expr.reshape(-1, 1)  # Reshape for HMM input (each cell's data is one sample)
        hidden_states = model.predict(cell_expr)  # Predict the hidden states (CNV states) for the cell
        predicted_cnv = np.array([cnv_levels[state] for state in hidden_states])  # Map states to CNV levels
        predicted_cnv_all.append(predicted_cnv)  # Store the predicted CNV levels for this cell
        cnv_states_all.append(hidden_states)  # Store the predicted hidden states for this cell

    # Convert the list of predicted CNV levels to a numpy array
    predicted_cnv_all = np.array(predicted_cnv_all)

    if plots:
        # Plot a heatmap showing the predicted CNV levels across cells and genes
        plt.figure(figsize=(12, 8))
        plt.imshow(predicted_cnv_all, aspect='auto', cmap='coolwarm', interpolation='none')
        plt.colorbar(label='CNV Level')
        plt.xlabel("Genes")
        plt.ylabel("Cells")
        plt.title(f"Predicted {cell_type} CNV Heatmap")
        plt.show()

        # Plot the distribution of normalized expression values to visualize the data spread
        plt.hist(expX.flatten(), bins=100)
        plt.title("Filtered + normalized expression distribution")
        plt.show()

    adata_filtered.layers['cnv_hmm'] = predicted_cnv_all

    # Add row wise means of hmm cnv matrix as cnv scores to adata_filtered.obs
    adata_filtered.obs['hmm_cnv_score'] = predicted_cnv_all.mean(axis=1)

    # Decide a threshold margin above and below hmm_cnv_score=1 to decide if it's a CNV or not
    threshold_margin = 0.05
    gain_mask = adata_filtered.obs['hmm_cnv_score']>(1+threshold_margin)
    loss_mask = adata_filtered.obs['hmm_cnv_score']<(1-threshold_margin)

    adata_filtered.obs['hmm_cnv'] = np.select([gain_mask,loss_mask], ['gain', 'loss'], default='')

    for cell_idx in range(predicted_cnv_all.shape[0]):
        cnv_type = adata_filtered.obs['hmm_cnv'][cell_idx]
        if cnv_type == 'gain':
            # Find gene indices that stay longest in 1.5
            value = 1.5
        elif cnv_type == 'loss':
            # Find gene indices that stay longest in 0.5
            value = 0.5
        else:
            continue

        # Keep finding longest region until start and end are on the same chromosome
        # In case the found indices lie on 2 different chromosomes
        chr_start,chr_end = 'a','b'
        start_idx = 0
        end_idx = len(predicted_cnv_all[cell_idx])

        while chr_start != chr_end:
            # Get indices of longest repeated value
            start_idx,end_idx = get_indices_repeated_value(predicted_cnv_all[cell_idx][start_idx:end_idx],value)

            # Get chromosome and region information for indices
            chr_start = adata_filtered.var['chromosome'][start_idx]
            chr_end = adata_filtered.var['chromosome'][end_idx]
            start = adata_filtered.var['start'][start_idx]
            end = adata_filtered.var['end'][end_idx]

        adata_filtered.obs['hmm_cnv'][cell_idx] = f'{chr_start}:{start}-{end} ({cnv_type})'

    return adata_filtered

def get_indices_repeated_value(arr,value):
    """ Finds the start and end index of the genomic region contained within a single chromosome that stays in a specified state the longest.

    Args:
        arr (list): list of HMM states per gene
        value (float): HMM state to target depending on gain (1.5) or loss (0.5)

    Returns:
        int tuple: start and end index
    """    
    max_len = 0
    start_index = -1
    end_index = -1
    index = 0

    for val, group in itertools.groupby(arr):
        group_list = list(group)
        length = len(group_list)
        if val == value and length > max_len:
            max_len = length
            start_index = index
            end_index = index + length - 1
        index += length

    return start_index, end_index

def extract_cnv_info(cnv_annotation):
    """ Splits a CNV annotation into chromosome, start, end, and type.

    Args:
        cnv_annotation (str): chr:start-end (type)

    Returns:
        str tuple: chromosome, start, end, type
    """    
    if cnv_annotation=='gain' or cnv_annotation=='loss':
        return [cnv_annotation]
    parts = cnv_annotation.split(" ")
    if len(parts) > 1:
        cnv_chr = parts[0].split(":")[0]
        cnv_start = parts[0].split(':')[1].split('-')[0]
        cnv_end = parts[0].split(':')[1].split('-')[1]
        cnv_type = parts[-1].strip("()")
        if cnv_type == '1' or cnv_type == '0':
            cnv_type = 'loss'
        if cnv_type == '4':
            cnv_type = 'gain'
        return cnv_chr,cnv_start,cnv_end,cnv_type
    return ''

def assess_predicted_cnvs(adata,prediction_annotation='hmm_cnv',truth_annotation='simulated_cnvs'):
    """ Prints precision, recall, accuracy, and F1 score of predicted CNV types

    Args:
        adata (AnnData): AnnData object with both simulated CNVs and predicted CNVs annotated
        prediction_annotation (str, optional): .obs column for predictions cnvs. Defaults to 'hmm_cnv'.
        truth_annotation (str, optional): .obs column for simulated cnvs. Defaults to 'simulated_cnvs'.
    """    
    truth = [extract_cnv_info(cnv)[-1] if cnv!='' else '' for cnv in adata.obs[truth_annotation]]
    prediction = [extract_cnv_info(cnv)[-1] if cnv!='' else '' for cnv in adata.obs[prediction_annotation]]
    df = pd.DataFrame({'truth':truth,'prediction':prediction})
    
    TP = ((df['truth'] == df['prediction']) & (df['truth'].isin(['gain', 'loss']))).sum()
    FP = ((df['truth'] == '') & (df['prediction'].isin(['gain', 'loss']))).sum()
    FN = ((df['truth'].isin(['gain', 'loss'])) & (df['prediction'] == '')).sum()
    TN = ((df['truth'] == '') & (df['prediction'] == '')).sum()

    precision = TP/(TP+FP)
    recall = TP/(TP+FN)
    accuracy = (TP+TN)/(TP+FN+TN+FP)
    F1_score = 2*precision*recall/(precision+recall)

    print(f'Precision: {precision}')
    print(f'Recall: {recall}')
    print(f'Accuracy: {accuracy}')
    print(f'F1 score: {F1_score}')