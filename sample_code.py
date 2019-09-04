#coding: utf-8

import numpy as np
import pandas as pd
import scEntropy.scEntropy as scEntropy
import matplotlib.pyplot as plt
import os

base_dir = os.getcwd()
# Get data dir
data_dir = os.path.join(base_dir, 'data')
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

# This function is used to process the raw HNSCC data, and then generate HNSCC_data_info.txt, HNSCC_data_content.txt,
# HNSCC_gene_expression.csv and HNSCC_cell_category.csv files.
# Since the raw file is already log2-transformed, we don't do this transformation again
def process_raw_data():

    # Read raw data
    raw_file_name = 'HNSCC_all_data.txt'
    raw_file_path = os.path.join(data_dir, raw_file_name)
    df_info = pd.read_csv(raw_file_path, sep='\t', quotechar='\'', dtype=str, nrows=5)
    df_content = pd.read_csv(raw_file_path, sep='\t', quotechar='\'', dtype=str).iloc[5:, :]

    # Output info file and content file
    # info_path = os.path.join(data_dir, 'HNSCC_data_info.txt')
    # df_info.to_csv(info_path, sep='\t', index=False, na_rep='')
    # content_path = os.path.join(data_dir, 'HNSCC_data_content.txt')
    # df_content.to_csv(content_path, sep='\t', index=False, na_rep='')

    # data type transformation
    df_content = df_content.iloc[:, 1:]
    df_content = df_content.applymap(lambda x: np.float(x))
    df_content = df_content.fillna(0)
    # The raw data has already log2 transformed, so we don't need to log2 transform again
    assert df_content.max().max() < 20

    # Output gene expression file
    expression_path = os.path.join(data_dir, 'HNSCC_gene_expression.csv')
    df_content.to_csv(expression_path, sep=',', index=False)

    # Generate cell category file(the cell categories in df_info have 3 types: normal, malignant, unknown)
    df_category = pd.DataFrame(columns=list(df_content.columns))
    # unknown type: -1, normal: 0, malignant: 1
    df_category.loc[0] = -1 * np.ones(df_category.shape[1], dtype=np.int64)
    cancer_class = (df_info.iloc[2:3, 1:] == '1').values[0]
    cancer_cell = df_content.loc[:, cancer_class].columns
    non_cancer_class = (df_info.iloc[3:4, 1:] == '1').values[0]
    normal_cell = df_content.loc[:, non_cancer_class].columns
    df_category.loc[0, cancer_cell] = 1
    df_category.loc[0, normal_cell] = 0

    # Output cell category file
    category_path = os.path.join(data_dir, 'HNSCC_cell_category.csv')
    df_category.to_csv(category_path, sep=',', index=False)


# This is a sample function to show how to use scEntropy with a predefined reference cell(RC)
def scEntropy_with_pre_ref():
    data_file_path = os.path.join(data_dir, 'HNSCC_gene_expression.csv')
    df_data = pd.read_csv(data_file_path)

    # Get scEntropy index array by the default settings
    scEntropy_arr = scEntropy.scEntropy(df_data, option='predefined')

    category_file_path = os.path.join(data_dir, 'HNSCC_cell_category.csv')
    df_category = pd.read_csv(category_file_path)
    cancer_class = (df_category.loc[0] == 1).values
    normal_class = (df_category.loc[0] == 0).values

    # Get scEntropy array for cancer cells and normal cells separately
    cancer_scEntropy_arr = scEntropy_arr[cancer_class]
    normal_scEntropy_arr = scEntropy_arr[normal_class]

    # The histogram boundaries and step length can be modified
    # We only generate histogram on [0, 5] range, using step=0.01.
    score_arr, _ = np.histogram(scEntropy_arr, bins=np.arange(0, 5.01, 0.01))
    cancer_score_arr, _ = np.histogram(cancer_scEntropy_arr, bins=np.arange(0, 5.01, 0.01))
    normal_score_arr, bin_edges = np.histogram(normal_scEntropy_arr, bins=np.arange(0, 5.01, 0.01))
    bin_edges = (bin_edges[1:] + bin_edges[:-1]) / 2

    # Visualization
    fig = plt.figure(figsize=(8, 8))
    plt.plot(bin_edges, score_arr, c='k')
    plt.plot(bin_edges, cancer_score_arr, c='r')
    plt.plot(bin_edges, normal_score_arr, c='b')
    plt.xlim([4, 5])
    # plt.ylim([0.45, 1])
    plt.xlabel('scEntropy')
    plt.ylabel('Density')
    fig.savefig('scEntropy_with_pre_ref.png', dpi=300)
    plt.show()
    return


# This is a sample function to show how to use scEntropy with RCSA. We call this scEGEMM framework.
# The scEGMM framework use RCSA to identify the intrinsic reference cell(IRC) in a data set, then the scEntropy index for each cell
# is calculated relative to the IRC.
def scEntropy_with_RCSA():
    data_file_path = os.path.join(data_dir, 'HNSCC_gene_expression.csv')
    df_data = pd.read_csv(data_file_path)

    # Get scEntropy index array by the default settings
    scEntropy_arr = scEntropy.scEntropy(df_data, option='RCSA')

    category_file_path = os.path.join(data_dir, 'HNSCC_cell_category.csv')
    df_category = pd.read_csv(category_file_path)
    cancer_class = (df_category.loc[0] == 1).values
    normal_class = (df_category.loc[0] == 0).values

    # Get scEntropy array for cancer cells and normal cells separately
    cancer_scEntropy_arr = scEntropy_arr[cancer_class]
    normal_scEntropy_arr = scEntropy_arr[normal_class]

    # The histogram boundaries and step length can be modified
    # We only generate histogram on [0, 5] range, using step=0.01.
    score_arr, _ = np.histogram(scEntropy_arr, bins=np.arange(0, 5.01, 0.01))
    cancer_score_arr, _ = np.histogram(cancer_scEntropy_arr, bins=np.arange(0, 5.01, 0.01))
    normal_score_arr, bin_edges = np.histogram(normal_scEntropy_arr, bins=np.arange(0, 5.01, 0.01))
    bin_edges = (bin_edges[1:] + bin_edges[:-1]) / 2

    # Visualization
    fig = plt.figure(figsize=(8, 8))
    plt.plot(bin_edges, score_arr, c='k')
    plt.plot(bin_edges, cancer_score_arr, c='r')
    plt.plot(bin_edges, normal_score_arr, c='b')
    plt.xlim([4, 5])
    # plt.ylim([0.45, 1])
    plt.xlabel('scEntropy')
    plt.ylabel('Density')
    fig.savefig('scEntropy_with_RCSA.png', dpi=300)
    plt.show()
    return


def main():
    process_raw_data()
    scEntropy_with_pre_ref()
    scEntropy_with_RCSA()
    return
