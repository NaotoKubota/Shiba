import argparse
import os
import sys
import pandas as pd
import scipy.stats as stats
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer


def get_args():
    ## Get arguments from command line

    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = "pca.py: Principal Component Analysis for matrix of gene expression and splicing"
    )

    parser.add_argument('--input-tpm', type=str, help='Input TPM file')
    parser.add_argument('--input-psi', type=str, help='Input PSI file')
    parser.add_argument('-g', '--genes', type=int, help='Number of highly-variable genes or splicing events to calculate PCs')
    parser.add_argument('-o', '--output', type=str, help='Output directory')

    args = parser.parse_args()

    return args


def load_tpm_table(tpm_file: str) -> pd.DataFrame:
    '''
    Load TPM table from input file

    Args:
    - tpm_file (str): input file containing TPM values

    Returns:
    - tpm_df (pd.DataFrame): dataframe containing TPM values
    '''

    tpm_df = pd.read_csv(tpm_file, sep="\t", index_col=0)
    # Drop rows with all zeros
    tpm_df = tpm_df.loc[(tpm_df != 0).any(axis=1)]
    # Drop rows with NaN
    tpm_df = tpm_df.dropna()

    return tpm_df


def load_psi_table(psi_file: str) -> pd.DataFrame:
    '''
    Load PSI table from input file

    Args:
    - psi_file (str): input file containing PSI values

    Returns:
    - psi_df (pd.DataFrame): dataframe containing PSI values
    '''

    psi_df = pd.read_csv(psi_file, sep="\t", index_col=0)
    psi_df = psi_df.drop(columns = ["pos_id"])

    # KNN imputation when psi_df has less than 6000 rows without NaN values
    if psi_df.dropna().shape[0] < 6000:
        print("Warning: PSI table has less than 6000 rows without NaN values. Performing KNN imputation...", file=sys.stderr)
        psi_df = psi_df.dropna(axis=0, thresh=psi_df.shape[1]*0.5) # Drop rows with more than 50% NaN
        print("Number of rows after dropping rows with more than 50% NaN: ", psi_df.shape[0], file=sys.stderr)
        # Check if there are any columns with all NaN
        column_all_nan = psi_df.columns[psi_df.isnull().all()]
        print("Columns with all NaN: ", column_all_nan, file=sys.stderr)
        print("Dropping columns with all NaN...", file=sys.stderr)
        psi_df = psi_df.drop(columns=column_all_nan)
        imputer = KNNImputer(n_neighbors=5)
        psi_df = pd.DataFrame(imputer.fit_transform(psi_df), index=psi_df.index, columns=psi_df.columns)
        print("KNN imputation done!", file=sys.stderr)
    else:
        # Drop rows with NaN
        psi_df = psi_df.dropna()

    return psi_df


def mtx2pca(df, genes) -> pd.DataFrame:
    '''
    Perform PCA on the input dataframe

    Args:
    - df (pd.DataFrame): input dataframe
    - genes (int): number of highly-variable genes to calculate PCs

    Returns:
    - feature_df (pd.DataFrame): dataframe containing principal components
    - contribution_df (pd.DataFrame): dataframe containing the contribution of each principal component
    '''

    # Keep rows of top n highly-variable genes
    if df.shape[0] > genes:
        df = df.loc[df.var(axis=1).sort_values(ascending=False).index[:genes]]

    # Z-score normalization across samples
    normalized_df = df.T.apply(stats.zscore, ddof = 1)

	# PCA
    pca = PCA()
    pca.fit(normalized_df)

    # Feature
    feature = pca.transform(normalized_df)
    feature_df = pd.DataFrame(feature, columns=["PC{}".format(x + 1) for x in range(len(feature))])
    feature_df.index = df.columns

    # Contribution
    contribution_df = pd.DataFrame(pca.explained_variance_ratio_, index=["PC{}".format(x + 1) for x in range(len(feature))])

    return(feature_df, contribution_df)


def main():

    # Get arguments
    args = get_args()

    # Load input files
    print("Loading input files...", file=sys.stderr)
    tpm_df = load_tpm_table(args.input_tpm)
    psi_df = load_psi_table(args.input_psi)

    # Perform PCA
    print("Performing PCA...", file=sys.stderr)
    tpm_feature_df, tpm_contribution_df = mtx2pca(tpm_df, args.genes)
    psi_feature_df, psi_contribution_df = mtx2pca(psi_df, args.genes)

    # Save output
    print("Saving output...", file=sys.stderr)
    os.makedirs(args.output, exist_ok=True)
    tpm_feature_df.to_csv(os.path.join(args.output, "tpm_pca.tsv"), sep="\t")
    tpm_contribution_df.to_csv(os.path.join(args.output, "tpm_contribution.tsv"), sep="\t", header=False)
    psi_feature_df.to_csv(os.path.join(args.output, "psi_pca.tsv"), sep="\t")
    psi_contribution_df.to_csv(os.path.join(args.output, "psi_contribution.tsv"), sep="\t", header=False)

    print("Done!", file=sys.stderr)

if __name__ == '__main__':

    main()
