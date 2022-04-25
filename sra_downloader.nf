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
  storeDir "${params.outdir}/raws/"
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
    path "*_fastqc.html", emit: fastqc_html_report
    path "*_fastqc.zip", emit: fastqc_zip_report
  script:
    """
    fastqc ${fastas}
    """
}

process fastp {
  storeDir "${params.outdir}/trimmed"
  container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.2--hb7a2d85_2"
  input:
    path fastas
  output:
    path "*_snipped.fastq", emit: trimmed_fastq
    path "*.html", emit: fastp_html_report
    path "*.json", emit: fastp_json_report
  script:
    if (fastas instanceof List){
      """
      fastp -i ${fastas[0]} -I ${fastas[1]} -o ${fastas[0].getSimpleName()}_snipped.fastq -O ${fastas[1].getSimpleName()}_snipped.fastq -h ${fastas[0].getSimpleName()}_fastp.html -j ${fastas[0].getSimpleName()}_fastp.json
      """
    }
    else {
      """
      fastp -i ${fastas} -o ${fastas.getSimpleName()}_snipped.fastq -h ${fastas[0].getSimpleName()}_fastp.html -j ${fastas[0].getSimpleName()}_fastp.json
      """
    }
}

process multiqc {
  storeDir "${params.outdir}/multiqc"
  container "https://depot.galaxyproject.org/singularity/multiqc:1.9--py_1"
  input:
    path inreps
  output:
    path "*"
  script:
    """
    multiqc .
    """
}

workflow {
  sraresult = prefetch(params.accession)
  dumpedfastas =  dumpprocess(sraresult)
  fastp_out = fastp(dumpedfastas)
  fastqc_out = fastqc(dumpedfastas.combine(fastp_out.trimmed_fastq))
  multiqc( fastp_out.fastp_json_report.concat(fastqc_out.fastqc_zip_report).collect() )
}
