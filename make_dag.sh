snakemake --dag > figures/dag.tmp
cat figures/dag.tmp | dot -Tpng > figures/dag.png
