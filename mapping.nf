nextflow.enable.dsl = 2

params.outdir = "results"
params.large_index = false
params.reference = null

if (!params.reference) {
    error "ERROR: Please specify reference genome"
}

process bowtie2_index {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/bowtie2%3A2.4.5--py39hd2f7db1_2"
  input:
    path ref
  output:
    path "*", emit: indices
  script:
    """
    bowtie2-build ${ref} ${ref.getBaseName()}_index
    """
}

process bowtie_map {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/bowtie2%3A2.4.5--py39hd2f7db1_2"
  input:
//    path ref
    path sample
    path indices
  output:
    path "*"
  script:
  """
  bowtie2 -x ${indices[0].getSimpleName()} -U ${sample} -S ${indices[0].getSimpleName()}.sam
  """
}

workflow {
//  refchannel = channel.fromPath("${params.reference}/*.fasta")
  refchannel = channel.fromPath("${params.reference}")
  samplechannel = channel.fromPath("${params.sample}/*.fastq")
  bowtie = bowtie2_index(refchannel)
  bowtie_map(samplechannel, bowtie.indices)

}
