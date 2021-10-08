# RESA: single cell Recurrently Expressed SNV Analysis
Copyright:

## Summary
RESA detects somatic mutation with high precision directly from scRNA-seq data. The key assumption of this method is that cancer cells evolve in a clonal manner and thus expressed somatic mutations have cross-cell recurrence, whereas the noise and artefacts are likely distributed randomly with small probability of recurrence across the cell population. RESA is composed of three main steps: initial variant calling, filtering and labeling, and modeling and refinement.


## RESA workflow
![image](https://user-images.githubusercontent.com/8051136/136513663-8e0f5a8f-29d2-44d2-a7a4-5bed334c3124.png)


## Installation
### Dependencies
Perl (version >=5)

R >= 3.4

Python >=3.7

Python Packages:
pandas>=1.3.0, rpy2 >=2.9.4, imblearn >= 0.8.0, argparse, sys, os, sklearn, numpy

### Install RESA

## Running RESA
1. Variant calling


2. RESA-filter


3. RESA-refine


