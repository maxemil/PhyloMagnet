
/*

*/
process alignContigs {
    input:
    file translated_contigs
    set file(reference_alignment) from EggNOGAlignments

    output:
    file '*.aln' into aligned_contigs

    publishDir 'queries', mode: 'copy'

    """
    echo $reference_alignment
    #mafft-fftnsi --adjustdirection --thread %s --addfragments ${translated_contigs} ${reference_alignment}
    """
}
