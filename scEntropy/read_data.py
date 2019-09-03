#coding: utf-8

import numpy as np
import pandas as pd

# Get default ref from the dataframe(a data structure in numpy structure) of a gene expression matrix
# Parameters:
# 1. file_path: the file_path of a gene expression file
# 2. log: whether the gene expression value need log2-transformation or not.
#    The common raw data value type is counts, TPM, FPKM and RPKM.
#    The default option is two apply log2-transformation on the dataframe
# 3. header: Whether the data file has header(column names) or not
#    !!! Pay attention that we assume the first column of the data file is gene name/symbol column
def read_data(file_path, log=True, header=True):
    if header == True:
        df_data = pd.read_csv(file_path, sep=',', dtype=str)
    else:
        df_data = pd.read_csv(file_path, sep=',', header=None, dtype=str)

    # We assume the first column of the dataframe is gene name/symbol column, so columns ranging from 2 to the last are data value columns
    df_content = df_data.iloc[:, 1:]
    df_content = df_content.applymap(np.float)
    df_content = df_content.fillna(0)
    if log == True:
        df_content = df_content.applymap(lambda x: np.log2(1.0 + x))
    return df_content


# Filter some genes which make less contribution to the subsequent scEntropy calculation
# Parameters:
# 1. df_content: the dataframe of gene expression matrix without the gene name column
def gene_selection_by_variance(df_content):
    df_expression_variance = df_content.apply(lambda x: x.var(), axis=1)
    high_variance_gene_idx = df_expression_variance[df_expression_variance > 1].index
    return df_content.loc[high_variance_gene_idx, :]