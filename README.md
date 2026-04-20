# RECUR ![Tested](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/OrthoFinder/RECUR/main/badge-status.json)

Finding Recurrent Substitutions from Multiple Sequence Alignments

## Introduction
![RECUR method workflow](./docs/images/RECUR_workflow.png)

<div align="center">
  Figure 1: The RECUR workflow
</div>

The required input is either a protein or codon multiple sequence alignment (in FASTA format) and a defined outgroup species or clade. The output of RECUR is a list of recurrent amino acid substitutions (file suffix: `.recur.tsv`), that have occurred in the inferred phylogeny. Outputs of intermediate steps, i.e. model selection, tree inference, ancestral state reconstruction and site substitution matrices, can be found in the `.recur` output directory.

Comprehensive installation and usage documentation is available on the [RECUR homepage](https://orthofinder.github.io/RECUR/)

## Citations

*Elizabeth H J Robbins, Yi Liu, Steven Kelly*, **RECUR: identifying recurrent amino acid substitutions from multiple sequence alignments**, *Molecular Biology and Evolution*, Volume 43, Issue 2, February 2026, msag036  [![DOI](https://img.shields.io/badge/DOI-10.1093%2Fmolbev%2Fmsag036-blue)](https://doi.org/10.1093/molbev/msag036)

## Credits and Acknowledgements

This is a software developed by the [Steven Kelly Lab](http://www.stevekellylab.com/).

## Contributing

If you find a bug :bug:, please open a [bug report](https://github.com/).
If you have an idea for an improvement or new feature, please open a [feature request]().

