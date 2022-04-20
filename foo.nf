nextflow.enable.dsl=2

process count_codon {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path infile
    val codon
  output:
    path "${infile}_${codon}.count"
    path "grepresult"
  script:
    """
    grep -oi ${codon} ${infile} > grepresult || :
    cat grepresult | wc - > ${infile}_${codon}.count
    """
}

workflow {
  inchannel = channel.fromPath(params.infile)
  count_word(inchannel, params.codon)
}
