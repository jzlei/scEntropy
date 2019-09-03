# scEntropy
We propose **scEntropy method** and **scEGMM framework**(see [[1]](https://www.biorxiv.org/content/10.1101/678557v1)) to measure the order of the cellular transcriptome profile from single cell RNA-seq data.


# Getting started
scEntropy method is a reference based method to obtain scEntropy index for each cell in a data set. We provide two different options to generate scEntropy. One is to use predefined reference cell expression(RC), and then calculate scEntropy index for each cell relative to RC. The other is to apply scEGMM framework to automatically identify reference cell expression in a data set and accordingly calculate scEntropy for each cell.

Here are the examples in which we use a data set which contains thousands of cells from 18 head and neck squamous cell carcinoma(HNSCC) patients. The HNSCC data set can be downloaded from [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322). The raw data is processed first to obtain clean gene expression data. The details of data pre-processing are shown in [**sample_code.py**](https://github.com/jzlei/scEntropy/blob/master/sample_code.py).

The examples of two different scEntropy options are given as follows:

**scEntropy with predefined reference cell**

`
import scEntropy
scEntropy.
`


**scEGMM framework**

`
import scEntropy
scEntropy.
`


# Final words

**Installation:**

`pip install git+https://github.com/jzlei/scEntropy.git`

**Requirements:**

- [Numpy](https://github.com/numpy/numpy)
- [Scipy](https://github.com/scipy/scipy)
- [Pandas](https://github.com/pandas-dev/pandas)
- [Scikit-Learn](https://github.com/scikit-learn/scikit-learn)

**Papers:**

- [[1] Single-cell entropy to quantify the cellular transcriptome from single-cell RNA-seq data](https://www.biorxiv.org/content/10.1101/678557v1)

