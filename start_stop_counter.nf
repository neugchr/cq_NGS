nextflow.enable.dsl=2

process count_st_codons{
  publishDir "/tmp/st_count", mode: 'copy', overwrite: true
  input:
    path fastafile
    val codon
  output:
    //path "${fastafile}.counts", emit: counts
    path "${codon}.count"
  script:
    """
    grep -oi ${codon} ${fastafile} | wc -l > ${codon}.count
    """
}


workflow {
  filechannel = channel.fromPath("/home/cq/NGS/cq-examples/data/sequences.fasta")
  wordchannel = channel.fromList(["ATG", "TAA", "TAG", "TGA"])
  count_st_codons(filechannel, wordchannel)
}
