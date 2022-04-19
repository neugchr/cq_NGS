nextflow.enable.dsl=2

process split_file {
  input:
    path fastafile
  output:
    path "${fastafile}.line*"
  script:
    """
    split -l 2 ${fastafile} -d ${fastafile}.line
    """
}

process count_codon{
//  publishDir "/tmp/st_count", mode: 'copy', overwrite: true
  publishDir "/tmp/fasta_count", mode: 'copy', overwrite: true
  input:
    path sequence
    val codon
  output:
    //path "${fastafile}.counts", emit: counts
    path "${codon}.count"
  script:
    """
    grep -oi ${codon} ${sequence} | wc -l > ${codon}.count
    """
}


workflow {
  splitted_seqs = split_file("/home/cq/NGS/cq-examples/data/sequences.fasta").flatten()
  wordchannel = channel.fromList(["ATG", "TAA", "TAG", "TGA"])
  codon = count_codon(splitted_seqs, wordchannel)
  wordchannel.view()
  splitted_seqs.view()
}
/*
workflow {
  filechannel = channel.fromPath("/home/cq/NGS/cq-examples/data/sequences.fasta")
  wordchannel = channel.fromList(["ATG", "TAA", "TAG", "TGA"])
  count_st_codons(filechannel, wordchannel)
}
*/
