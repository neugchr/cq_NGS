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
  publishDir "/tmp/fasta_count", mode: 'copy', overwrite: true
  input:
    tuple path(sequence), val(codon)
  output:
    path "${codon}.count"
  script:
    """
    grep -oi ${codon} ${sequence} | wc -l > ${codon}.count
    """
}


workflow {
  splitted_seqs = split_file("/home/cq/NGS/cq-examples/data/sequences.fasta").flatten()
  wordchannel = channel.fromList(["ATG", "TAA", "TAG", "TGA"])
  combined_channel = splitted_seqs.combine(wordchannel)
  codon = count_codon(combined_channel)
}
