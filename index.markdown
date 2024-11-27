---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: default
title: Home
---

## Introduction
![RECUR method workflow](./assets/images/RECUR_workflow_figure.png)

<div align="center">
  Figure 1: The RECUR workflow
</div>

The required input is either a protein or codon multiple sequence alignment (in FASTA format) and a defined outgroup species or clade. The output of RECUR is a list of recurrent amino acid substitutions, that have occurred in the inferred phylogeny (file suffix: `.recur.tsv`). Outputs of intermediate steps, i.e. model selection, tree inference, ancestral state reconstruction and site substitution matrices, can be found in the .recur output directory.
