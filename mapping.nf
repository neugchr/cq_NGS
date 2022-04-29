nextflow.enable.dsl = 2

params.outdir = "results"
params.large_index = false
params.reference = "/home/cq/Documents/Module_3_NGS/cq_NGS/ON074439"

// bowtie 2 indexing fasta file ON074439
// bowtie 2 mapping reads onto index

process bowtie2_index {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/bowtie2%3A2.4.5--py39hd2f7db1_2"
  input:
    path ref
  output:
    path "*"
  script:
    """
    bowtie2-build ${ref} index
    """
}

workflow {
//  fastqchannel = channel.fromPath("${params.indir}/*.fastq").collect()
  refchannel = channel.fromPath("${params.reference}/*.fasta")
  bowtie2_index(refchannel)
}
