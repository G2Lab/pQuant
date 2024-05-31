snakemake \
    --configfile=config/config.yaml \
    --use-conda \
    --profile=simple/ \
    --keep-going \
    -s Snakefile \
    "$@"
