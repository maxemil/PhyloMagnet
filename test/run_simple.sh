function run_simple()
{
  nextflow run main.nf \
            -with-singularity PhyloMagnet.sif \
            --is_runs true \
            --fastq "test/*rpoB.fastq.gz" \
            --reference_packages "test/rpkgs/*" \
            --lineage "order" \
            --megan_eggnog_map eggnog.map \
            --cpus 2 \
            --is_runs true \
            --queries_dir test/queries \
            --reference_dir test/references \
            --phylo_method 'fasttree' \
            --align_method 'mafft-fftnsi' \
            -w test/work -resume

if $(grep -q "Enterobacterales\sTrue" test/queries/tree_decisions.txt)
then
  return 0
else
  return 1
fi
}

singularity pull --name PhyloMagnet.sif shub://maxemil/PhyloMagnet

run_simple
