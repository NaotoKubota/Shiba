# Installation

You don't have to build a local environment for running Shiba! We provide a Docker image of **Shiba**, which includes all the dependencies and the **Shiba** software itself. You can run them using [Docker](https://docs.docker.com/get-docker/)/[Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) or [Snakemake](https://snakemake.readthedocs.io/en/stable/).

## Docker

``` bash
# Pull the latest image
docker pull naotokubota/shiba:latest

# Login to the container
docker run -it --rm naotokubota/shiba:latest bash

# Run Shiba, for example, to see the help message
docker run -it --rm naotokubota/shiba:latest Shiba -h
```

## Singularity

``` bash
# Pull the latest image
singularity pull docker://naotokubota/shiba:latest

# Login to the container
singularity shell shiba_latest.sif

# Run Shiba, for example, to see the help message
singularity exec shiba_latest.sif Shiba -h
```

## Snakemake

You need to install [Snakemake](https://snakemake.readthedocs.io/en/stable/) and clone the **Shiba** GitHub repository on your system:

``` bash
# Clone the Shiba repository
git clone https://github.com/NaotoKubota/Shiba.git
```

And please make sure you have [Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html) installed on your system as the snakemake workflow uses Singularity to run each step of the pipeline.
