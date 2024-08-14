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
