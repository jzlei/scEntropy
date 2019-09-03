#coding: utf-8

from setuptools import setup, find_packages

install_requires = [
    'numpy>=1.12.1',
    'pandas>=0.19.0',
    'scikit-learn>=0.19.2',
    'scipy>=1.1.0'
]

setup(
    name = 'scEntropy',
    version = '1.0.0',
    keywords='bioinformatics, data processing, entropy, single-cell RNA sequencing data',
    description = 'This is a library for scRNA-seq data analysis by scEntropy',
    license = 'GNU GPL v3',
    url = 'https://github.com/jzlei/scEntropy',
    author = 'Jingxin Liu, You Song, Jinzhi Lei',
    author_email = 'jzlei@tsinghua.edu.cn',
    packages = find_packages(),
    include_package_data = True,
    platforms = 'any',
    install_requires = install_requires,
)