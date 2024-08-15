# SnakeScShiba Usage

``` bash
snakemake -s SnakeScShiba \
--configfile config.yaml \
--cores 32 \
--use-singularity \
--singularity-args "--bind $HOME:$HOME" \
--rerun-incomplete
```

Please check the [Manual](../manual/diff_splicing_sc.md/#1-prepare-inputs_1){ data-preview } to learn how to prepare the `config.yaml`.

<figure markdown="span">
	![SnakeScShiba rulegraph](https://github.com/NaotoKubota/Shiba/blob/mkdocs/img/SnakeScShiba_rulegraph.svg?raw=true){ width="500" align="center" }
	<figcaption>SnakeScShiba rulegraph</figcaption>
</figure>

This rulegraph was made by [snakevision](https://github.com/OpenOmics/snakevision).
