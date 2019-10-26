#coding: utf-8

import numpy as np
import pandas as pd
import scipy.stats as st
from sklearn.mixture import GaussianMixture
import time

# Get default ref from the dataframe(a data structure in numpy structure) of a gene expression matrix
# Parameters:
# 1. df_content: the dataframe of a gene expression matrix
def get_default_ref(df_content):
    ref_vec = df_content.apply(np.mean, axis=1)
    return ref_vec

# Calculate Shanno Entropy for the series(a data structure in numpy package) of a gene expression vector
# Parameters:
# 1. data_series: the series of a gene expression vector
# 2. left: the left boundary value of relative gene expression
# 3. right: the right boundary value of relative gene expression
# 4. step: the range interval of the histogram of data_series
def calc_entropy(data_series, left=-100, right=100, step=0.01):
    bins_interval = np.arange(left - step / 2, right + 3 * step / 2, step)
    bins_cnt, bin_edges = np.histogram(data_series.values, bins=bins_interval)
    bins_cnt = bins_cnt[bins_cnt > 0]
    prob_vec = bins_cnt / bins_cnt.sum()
    res_sum = (-prob_vec * np.log2(prob_vec)).sum()
    return res_sum

# RCSA
# We use RCSA to find the intrinsic reference cell in a data set
# Parameters:
# 1. df_content: the dataframe of a gene expression matrix
# 2. components: the number of gaussian components in scEGMM framework, we use 2 gaussian components as default
def RCSA(df_content, components=2):

    idx = 0
    ref_vec = get_default_ref(df_content)
    #Calculate the initial error
    error = ref_vec.apply(lambda x: x**2).sum() / ref_vec.shape[0]
    clf = GaussianMixture(n_components=components, covariance_type='full')

    while error > 1e-6:
        start = time.time()
        #Calculate the entropy vector
        series_entropy = df_content.sub(ref_vec, axis=0).apply(calc_entropy, axis=0)
        df_entropy = pd.DataFrame({'entropy': series_entropy})
        end = time.time()
        print('calc entropy finish! time={0}s'.format(end - start))
        #GMM fitting
        clf.fit(df_entropy)
        if (clf.means_[0, 0] < clf.means_[1, 0]):
            comp = 0
        else:
            comp = 1
        #Get GMM model coefficients
        gau_wt = clf.weights_[comp]
        gau_mean = clf.means_[comp]
        gau_sigma = np.sqrt(clf.covariances_[comp, 0])

        #Calculate the cell weight
        entropy_den = gau_wt * st.norm.pdf(series_entropy, loc=gau_mean, scale=gau_sigma)
        cell_weight = entropy_den / entropy_den.sum()

        def calc_ref(series_vec, weight=None):
            return (series_vec.values * weight).sum()

        pre_ref = ref_vec.copy()
        ref_vec = df_content.apply(calc_ref, axis=1, weight=cell_weight)
        error = pre_ref.sub(ref_vec).apply(lambda x: x ** 2).sum() / ref_vec.shape[0]
        idx += 1
        end2 = time.time()
        print('idx = ', idx, ', error=', error, 'time={0}s'.format(end2 - end))
    return ref_vec

# Generate scEntropy in two different ways
# We provide two options
# 1. predefined
# You need to calculate the reference cell expression level by yourself, the default predefined reference cell expression
# is the average expression of all cells in a data set. You could also use another reference cell expression level vector as input.
# 2. RCSA
# We use RCSA in our paper to obtain the intrinsic reference cell expression level vector
# After the above reference cell calculation process, the scEntropy relative to reference cell is then generated.
# ----------------------------------
# Parameters:
# 1. df_content: the dataframe of a gene expression matrix
# 2. ref_vec: the reference cell expression vector for df_content
def scEntropy(df_content, ref_vec=None, option='predefined'):
    if option == 'predefined':
        if ref_vec is None:
            ref_vec = get_default_ref(df_content)
    elif option == 'RCSA':
        ref_vec = RCSA(df_content)
    else:
        raise Exception('parameter error in scEntropy function')

    if isinstance(ref_vec, list) or isinstance(ref_vec, np.ndarray):
        ref_vec = pd.Series(ref_vec)

    if ref_vec.shape[0] != df_content.shape[0]:
        raise Exception('The dimension of ref_vec is not consistent with df_content!')

    series_entropy = df_content.sub(ref_vec, axis=0).apply(calc_entropy, axis=0)
    return series_entropy.values
