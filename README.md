# Benchmarking microbial profilers on the CAMI data

## Running the pipeline
The pipeline has been implemented in Snakemake and to run using the module system within the HPC [Computerome](https://www.computerome.dk/).

To start the pipeline using the [Snakefile](Snakefile):
```{bash}
snakemake --cores 40 --use-envmodules --profile profile/ --configfile profile/config.yaml
```

## Converting output of profilers into the CAMI format
To convert the results of mOTUs, MetaPhlAn, and bracken we used the tocami.py script from https://github.com/hzi-bifo/cami2_pipelines/blob/master/bin/tocami.py. 

Using the .mapstat files produced by KMA, we converted them into the CAMI format using the [`kma2cami.py` script](kma2cami.py).

## Visualizing the results
The binary classification metrics are visualized in the jupyter notebook [`Comparison.ipynb`](Comparison.ipynb)

## Data availability
The in-silico data was retrieved from the following urls:
* https://frl.publisso.de/data/frl:6425521/
* https://frl.publisso.de/data/frl:6425518/

CAMI formatted results are available on [Zenodo](https://zenodo.org/record/7923775).