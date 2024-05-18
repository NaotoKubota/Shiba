[![DOI](https://zenodo.org/badge/801784656.svg)](https://zenodo.org/doi/10.5281/zenodo.11214807)

# Shiba

An exon-centric computational pipeline for robust identification of differential RNA splicing.

<img src="img/Shiba_icon.png" width=50%>

## Overview

<img src="img/Shiba_overview.png" width=75%>

## Usage

***Shiba***

```bash
Shiba -i experiment.tsv -c config.txt
```

***SnakeShiba***, Snakemake-based workflow of Shiba

```bash
snakemake -s SnakeShiba --configfile config.yaml --cores 32 --use-singularity
```

***scShiba***, a single-cell RNA-seq version of Shiba

```bash
scShiba -i experiment.tsv -c config.txt
```

***SnakeScShiba***, Snakemake-based workflow of scShiba

```bash
snakemake -s SnakeScShiba --configfile config.yaml --cores 32 --use-singularity
```

See [the manual document](https://github.com/NaotoKubota/Shiba/blob/main/doc/MANUAL.md) for details.
