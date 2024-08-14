# Differential RNA splicing analysis with bulk RNA-seq data

!!! note

	You need to install a Docker image of **Shiba** (and clone the **Shiba** GitHub repository to run **SnakeShiba**). If you don't have them installed, please follow the instructions in the [Installation](../installation.md) section.

---

## Before you start

- Please perform mapping of RNA-seq reads to the reference genome and generate bam files by software such as [STAR](https://github.com/alexdobin/STAR) and [HISAT2](https://daehwankimlab.github.io/hisat2/).
- Please download a gene annotataion file of your interest in GTF format.

Download a mouse gene annotation file (Ensembl 102):

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

`config.txt`: A text file of the configuration.

``` text
##############################
### Config file for Shiba ####
##############################

## General

# Number of processors to use
NUM_PROCESS=16
# Reference GTF file
GTF=/path/to/workdir/Mus_musculus.GRCm38.102.gtf
# Output directory
OUTPUT=/path/to/workdir/output
# Detect unannotated splicing events
UNANNOTATED=true

## Step3: bam2junc.sh

# Minimum anchor length for extracting exon-exon junction reads
MINIMUM_ANCHOR_LENGTH=6
# Minimum intron size for extracting exon-exon junction reads
MINIMUM_INTRON_SIZE=70
# Minimum intron size for extracting exon-exon junction reads
MAXIMUM_INTRON_SIZE=500000
# Strand specificity of RNA library preparation for extracting exon-exon junction reads
STRAND=XS

## Step4: psi.py

# Just calculate PSI for each sample, not perform statistical tests
ONLY_PSI=false
# Just calculate PSI for each group, not perform statistical tests
ONLY_PSI_GROUP=false
# FDR for DSE detection
FDR=0.05
# Miminum PSI change for DSE detection
DELTA_PSI=0.1
# Reference group for DSE detection
REFERENCE_GROUP=Ref
# Alternative group for DSE detection
ALTERNATIVE_GROUP=Alt
# Minumum value of total reads for each junction
MINIMUM_READS=10
# Print PSI for individual samples to output files
INDIVIDUAL_PSI=true
# Perform Welch's t-test between reference and Alternative group
TTEST=true

## Skip steps

SKIP_STEP1=false # bam2gtf.sh
SKIP_STEP2=false # gtf2event.py
SKIP_STEP3=false # bam2junc.sh
SKIP_STEP4=false # psi.py
SKIP_STEP5=false # expression.sh
SKIP_STEP6=false # pca.py
SKIP_STEP7=false # plots.py
```

You can skip some steps by setting `SKIP_STEP*` to `true`.

### 2. Run

Docker:

``` bash
cp experiment.tsv config.txt /path/to/workdir
cd /path/to/workdir
docker run --rm -v $(pwd):$(pwd) naotokubota/shiba:latest \
Shiba -i /path/to/workdir/experiment.tsv -c /path/to/workdir/config.txt
```

Singularity:

``` bash
cp experiment.tsv config.txt /path/to/workdir
singularity exec shiba_latest.sif \
Shiba -i /path/to/workdir/experiment.tsv -c /path/to/workdir/config.txt
```

!!! note

	When you use Singularity, you do not need to bind any paths as it automatically binds some paths in the host system to the container. In the default configuration, the system default bind points are `$HOME`, `/sys:/sys`, `/proc:/proc`, `/tmp:/tmp`, `/var/tmp:/var/tmp`, `/etc/resolv.conf:/etc/resolv.conf`, `/etc/passwd:/etc/passwd`, and `$PWD`. If files needed to be accessed are not in these paths, you can use the `--bind` option to bind the files to the container.

---

## SnakeShiba

A snakemake-based workflow of **Shiba**. This is useful for running **Shiba** on a cluster. Snakemake automatically parallelizes the jobs and manages the dependencies between them.

### 1. Prepare inputs

`experiment.tsv`: A tab-separated text file of sample ID, path to fastq files, and groups for differential analysis. This is the same as the input for **Shiba**.

`config.yaml`: A yaml file of the configuration.

``` yaml
workdir:
  /path/to/workdir
container:
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
only_psi:
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
  False
```

You can generate a file of splicing analysis results in excel format by setting `excel` to `True`.

### 2. Run

Please make sure that you have installed Snakemake and Singularity and cloned the **Shiba** repository on your system.

``` bash
snakemake -s /path/to/Shiba/SnakeShiba \
--configfile config.yaml \
--cores 16 \
--use-singularity \
--rerun-incomplete
```
