# scEntropy
We propose **scEntropy method** and **scEGMM framework**(see [[1]](https://www.biorxiv.org/content/10.1101/678557v1)) to measure the order of the cellular transcriptome profile from single cell RNA-seq data.


# Getting started
scEntropy method is a reference based method to obtain scEntropy index for each cell in a data set. We provide two different options to generate scEntropy. One is to **use predefined reference cell expression(RC)**, and then calculate scEntropy index for each cell relative to RC. The other is to **apply scEGMM framework to automatically identify reference cell** expression in a data set and accordingly calculate scEntropy for each cell.

Here are the examples in which we use a data set which contains thousands of cells from 18 head and neck squamous cell carcinoma(HNSCC) patients. The HNSCC data set can be downloaded from [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322). The raw data is processed first to obtain clean gene expression data. The details of data pre-processing are shown in [**sample_code.py**](https://github.com/jzlei/scEntropy/blob/master/sample_code.py).

The examples of two different scEntropy options are given as follows:

**scEntropy with predefined reference**

```python
# You need to use `pip install git+https://github.com/jzlei/scEntropy.git` to install this package
import scEntropy.scEntropy as scEntropy 
import pandas as pd
import os
data_file_path = os.path.join('data', 'HNSCC_gene_expression.csv')
df_data = pd.read_csv(data_file_path)
scEntropy.scEntropy(df_data, ref_vec=None, option='predefined')
```

If you don't pass any `ref_vec` as reference cell(RC), the predefined procedure will use the average gene expression of all cells to replace the `ref_vec`.
But if you pass a `ref_vec` which is not `None` type, the predefined procedure will use this passed reference cell to generate scEntropy index.

**scEGMM framework**

```python
# You need to use `pip install git+https://github.com/jzlei/scEntropy.git` to install this package
import scEntropy.scEntropy as scEntropy
import pandas as pd
import os
data_file_path = os.path.join('data', 'HNSCC_gene_expression.csv')
df_data = pd.read_csv(data_file_path)
scEntropy.scEntropy(df_data, option='RCSA')
```

If you use scEGMM framework, you don't need to pass any reference cell, since the scEGMM framework will automatically impose the RCSA to identify the intrinsic reference cell in a data set. And then scEntropy index will be generated.


To see the details of the above two different samples, please check the function `scEntropy_with_pre_ref` and `scEntropy_with_RCSA` in [**sample_code.py**](https://github.com/jzlei/scEntropy/blob/master/sample_code.py).

# Final words

Here are some useful information to install the scEntropy package and check details of scEntropy method.

**Installation:**

`pip install git+https://github.com/jzlei/scEntropy.git`

**Requirements:**

- [Numpy](https://github.com/numpy/numpy)
- [Scipy](https://github.com/scipy/scipy)
- [Pandas](https://github.com/pandas-dev/pandas)
- [Scikit-Learn](https://github.com/scikit-learn/scikit-learn)

**Papers:**

- [[1] Single-cell entropy to quantify the cellular transcriptome from single-cell RNA-seq data](https://www.biorxiv.org/content/10.1101/678557v1)

