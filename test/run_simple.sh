function run_simple()
{
  nextflow run main.nf \
            -with-singularity PhyloMagnet.simg \
            --reference_classes test/EggNOG.txt \
            --is_runs true \
            --fastq test/rpoB.fastq.gz \
            --local_ref test/rpoB_translation.fasta \
            --lineage "Enterobacterales" \
            --taxonomy_level_trees "order" \
            --megan_eggnog_map eggnog.map \
            --cpus 2 \
            --is_runs true \
            --queries_dir test/queries \
            --reference_dir test/references \
            -w test/work -resume

if $(grep -q "Enterobacterales\sTrue" test/queries/tree_decisions.txt)
then
  return 0
else
  return 1
fi
}

singularity pull --name PhyloMagnet.simg shub://maxemil/PhyloMagnet:latest

run_simple
