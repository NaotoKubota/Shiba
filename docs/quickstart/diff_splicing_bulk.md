# Differential RNA splicing analysis with bulk RNA-seq data

!!! note

	You need to install a Docker image of **Shiba** (and clone the **Shiba** GitHub repository to run **SnakeShiba**). If you don't have them installed, please follow the instructions in the [Installation](../installation.md) section.

---

## Before you start

- Perform mapping of RNA-seq reads to the reference genome and generate bam files by software such as [STAR](https://github.com/alexdobin/STAR) and [HISAT2](https://daehwankimlab.github.io/hisat2/).
- Download a gene annotataion file of your interest in GTF format.

Here is an example code for downloading a mouse gene annotation file (Ensembl 102):

``` bash
wget https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz
gzip -d Mus_musculus.GRCm38.102.gtf.gz
```

---

## Shiba

### 1. Prepare inputs

`experiment.tsv`: A tab-separated text file of sample ID, path to bam files, and groups for differential analysis.

``` text
sample<tab>bam<tab>group
sample_1<tab>/path/to/workdir/bam/sample_1.bam<tab>Ref
sample_2<tab>/path/to/workdir/bam/sample_2.bam<tab>Ref
sample_3<tab>/path/to/workdir/bam/sample_3.bam<tab>Ref
sample_4<tab>/path/to/workdir/bam/sample_4.bam<tab>Alt
sample_5<tab>/path/to/workdir/bam/sample_5.bam<tab>Alt
sample_6<tab>/path/to/workdir/bam/sample_6.bam<tab>Alt
```

Please put bam files with their index files (`.bai`) in the `path/to/workdir/bam` directory and replace `<tab>` with a tab character.

`config.yaml`: A yaml file of the configuration.

``` yaml
workdir:
  /path/to/workdir
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf
experiment_table:
  /path/to/experiment.tsv
unannotated:
  True

# Junction read filtering
minimum_anchor_length:
  6
minimum_intron_length:
  70
maximum_intron_length:
  500000
strand:
  XS

# PSI calculation
only_psi:
  False
only_psi_group:
  False
fdr:
  0.05
delta_psi:
  0.1
reference_group:
  Ref
alternative_group:
  Alt
minimum_reads:
  10
individual_psi:
  True
ttest:
  True
excel:
  True
```

You can generate a file of splicing analysis results in excel format by setting `excel` to `True`.

### 2. Run

Docker:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
cd /path/to/workdir
docker run --rm -v $(pwd):$(pwd) naotokubota/shiba:latest \
python /opt_shiba/Shiba/shiba.py -p 32 /path/to/workdir/config.yaml
```

Singularity:

``` bash
cp experiment.tsv config.yaml /path/to/workdir
singularity exec shiba_latest.sif \
python /opt_shiba/Shiba/shiba.py -p 32 /path/to/workdir/config.yaml
```

!!! note

	When you use Singularity, you do not need to bind any paths as it automatically binds some paths in the host system to the container. In the default configuration, the system default bind points are `$HOME`, `/sys:/sys`, `/proc:/proc`, `/tmp:/tmp`, `/var/tmp:/var/tmp`, `/etc/resolv.conf:/etc/resolv.conf`, `/etc/passwd:/etc/passwd`, and `$PWD`. If files needed to be accessed are not in these paths, you can use the `--bind` option to bind the files to the container.

---

## SnakeShiba

A snakemake-based workflow of **Shiba**. This is useful for running **Shiba** on a cluster. Snakemake automatically parallelizes the jobs and manages the dependencies between them.

### 1. Prepare inputs

`experiment.tsv`: A tab-separated text file of sample ID, path to fastq files, and groups for differential analysis. This is the same as the input for **Shiba**.

`config.yaml`: A yaml file of the configuration. This is the same as the configuration for **Shiba** but with the addition of the `container` field and without the `only_psi` and `only_psi_group` fields as they are not supported in **SnakeShiba**.

``` yaml
workdir:
  /path/to/workdir
container: # This field is required for SnakeShiba
  docker://naotokubota/shiba
gtf:
  /path/to/Mus_musculus.GRCm38.102.gtf
experiment_table:
  /path/to/experiment.tsv
unannotated:
  True

# Junction read filtering
minimum_anchor_length:
  6
minimum_intron_length:
  70
maximum_intron_length:
  500000
strand:
  XS

# PSI calculation
fdr:
  0.05
delta_psi:
  0.1
reference_group:
  Ref
alternative_group:
  Alt
minimum_reads:
  10
individual_psi:
  True
ttest:
  True
excel:
  True
```

You can generate a file of splicing analysis results in excel format by setting `excel` to `True`.

### 2. Run

Please make sure that you have installed Snakemake and Singularity and cloned the **Shiba** repository on your system.

``` bash
snakemake -s /path/to/Shiba/snakeshiba.smk \
--configfile config.yaml \
--cores 16 \
--use-singularity \
--rerun-incomplete
```
