# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 21:03:46 2021

@author: janes
"""

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="infect-net-inference",
    version="0.0.2",
    author="Min Xu and Hwai-Liang Tung",
    author_email="hwai-liang_tung@brown.edu",
    url="https://github.com/hltung/infection/tree/adjust-k",
    description="Simulate an infection process on a graph as well as create a credible set for the patient zero given a set of infected nodes on a graph.",
    py_modules=["infect_tools_noisy_conv", "tree_tools"],
    package_dir={'':'src'},
    classifiers=[],
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires = [""],
    extras_require = {
        "dev": [
            "pytest>=3.7",
            ],
        },
)