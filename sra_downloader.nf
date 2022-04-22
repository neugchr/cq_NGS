nextflow.enable.dsl = 2

//  container "https://depot.galaxyproject.org/singularity/"

process prefetch {
  storeDir "${params.outdir}"
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"
  input:
    val accession
  output:
    path "${accession}/${accession}.sra"
  script:
    """
    prefetch $accession
    """
}

process dumpprocess {
  storeDir "${params.outdir}"
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"
  input:
    path sraresult
  output:
    path "${sraresult.getBaseName()}*.fastq"   // could also be getSimpleName()
  script:
    """
    fasterq-dump.2.11.0 ${sraresult}
    """
}

process fastqc {
//  publishDir "${params.outdir}", mode: "copy", overwrite: true
  storeDir "${params.outdir}"
  container "https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1"
  input:
    path fastas
  output:
    path "*_fastqc.*"
  script:
    """
    fastqc ${fastas}
    """
}

workflow {
  sraresult = prefetch(params.accession)
  dumpedfastas = dumpprocess(sraresult)
  fastqc(dumpedfastas)
}
