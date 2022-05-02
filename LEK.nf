nextflow.enable.dsl = 2

params.indir = "/home/cq/Documents/Module_3_NGS/cq_NGS/LEK/"
params.outdir = "${params.indir}/results"
params.resistance_db = "/home/cq/Documents/Module_3_NGS/cq_NGS/LEK/CARD_v3.0.8_SRST2.fasta"
params.only_resistance = false


process fastp {
  storeDir "${params.outdir}/trimmed"
  container "https://depot.galaxyproject.org/singularity/fastp%3A0.23.2--hb7a2d85_2"
  input:
    path fastas
  output:
    path "${fastas.getSimpleName()}_snipped.fastq", emit: trimmed_fastq
    path "${fastas.getSimpleName()}_fastp.html", emit: fastp_html_report
    path "${fastas.getSimpleName()}_fastp.json", emit: fastp_json_report
  script:
    """
    fastp -i ${fastas} -o ${fastas.getSimpleName()}_snipped.fastq -h ${fastas.getSimpleName()}_fastp.html -j ${fastas.getSimpleName()}_fastp.json
    """
}

process fastqc {
  storeDir "${params.outdir}/fastqc"
  container "https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1"
  input:
    path fastas
  output:
    path "${fastas.getSimpleName()}_fastqc.html", emit: fastqc_html_report
    path "${fastas.getSimpleName()}_fastqc.zip", emit: fastqc_zip_report
  script:
    """
    fastqc ${fastas}
    """
}

process srst2 {
  storeDir "${params.outdir}/srst"
  container "https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2"
  input:
    path sample
  output:
    path "${sample.getBaseName()}*", emit: resistances
    path "${sample.getBaseName()}*__genes**", emit: genes
    path "${sample.getBaseName()}*__fullgenes*", emit: fullgenes
  script:
    """
    srst2 --input_se ${sample} --output ${sample.getBaseName()} --log --gene_db ${params.resistance_db}
    """
}

process aggregate {
  storeDir "${params.outdir}/aggregate_report"
  container "https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2"
  input:
    path genes
  output:
    path "*"
  script:
    """
    srst2 --prev_output ${genes} --output genes.report
    """
}

workflow {
  input = channel.fromPath("${params.indir}/*.fastq").collect().flatten()
  fastp_out = fastp(input)
  fastqc_out = fastqc(fastp_out.trimmed_fastq)
  srst2_input = fastp_out.trimmed_fastq
  srst2_out = srst2(srst2_input.flatten())
  // srst2_out.genes.collect().view()
  aggregate( srst2_out.genes.collect())
}
