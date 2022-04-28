nextflow.enable.dsl = 2

params.hashlen = 51
params.outdir = "results"

params.quastref = false

process velvet{
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/velvet:1.2.10--h7132678_5"
  input:
    path fastq
    val hashlen
  output:
    path "velvetdir"
    path "velvetdir/velvetcontigs.fa", emit: velvetcontigs
  script:
    if(fastq instanceof List){
      """
      velveth velvetdir ${hashlen} -shortPaired -fastq -separate ${fastq}
      velvetg velvetdir
      mv velvetdir/contigs.fa velvetdir/velvetcontigs.fa
      """
    } else {
      """
      velveth velvetdir ${hashlen} -fastq -short ${fastq}
      velvetg velvetdir
      mv velvetdir/contigs.fa velvetdir/velvetcontigs.fa
      """
    }
}

process quast{
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"
  input:
    path infile
  output:
    path "quast_results"
  script:
    """
    quast ${infile}
    """
}

process quast_ref{
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2"
  input:
    path infile
    path ref
  output:
    path "quast_results"
  script:
    println infile
    println ref
    """
    quast ${infile} -r ${ref}
    """
}

process spades {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/spades%3A3.15.4--h95f258a_0"
  input:
    path fastq
  output:
    path "spades_out"
    path "spades_out/spadescontigs.fa", emit: spadescontigs
  script:
  if(fastq instanceof List){
    """
    spades.py --only-assembler -1 ${fastq[0]} -2 ${fastq[1]} -o spades_out
    mv spades_out/contigs.fasta spades_out/spadescontigs.fa
    """
  } else {
    """
    spades.py --only-assembler -s ${fastq} -o spades_out
    mv spades_out/contigs.fasta spades_out/spadescontigs.fa
    """
  }
}

workflow {
  fastqchannel = channel.fromPath("${params.indir}/*.fastq").collect()

  vout = velvet(fastqchannel, params.hashlen)
  spds_out = spades(fastqchannel)
  all_contigs = vout.velvetcontigs.concat(spds_out.spadescontigs).collect()

  if(params.quastref){
    reference = channel.fromPath(params.quastref)
    quast_ref(all_contigs, reference)
  } else {
    quast(all_contigs)
  }

}
