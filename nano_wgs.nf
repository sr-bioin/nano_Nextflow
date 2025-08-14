#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Bacterial Whole Genome Sequencing Analysis Pipeline
 * FastQC → NanoPlot → MultiQC → Flye → QUAST → Prokka → Abricate
 * Only Flye uses an external conda environment
 */

params.reads   = "data/*.fastq.gz"
params.outdir  = "results"

// params.flye_envpath = "/home/user/anaconda3/envs/flye_env"     // external Flye env


// Channels
Channel.fromPath(params.reads)
       .map { file -> tuple(file.baseName.replaceAll(/\.fastq.*/, ""), file) }
       .set { reads_ch }

// ==========================
// 1. FastQC
// ==========================
process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.{html,zip}"), emit: reports

   // conda 'bioconda::fastqc=0.12.1'

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}

// ==========================
// 2. NanoPlot
// ==========================
process nanoplot {
    tag "$sample_id"
    publishDir "${params.outdir}/nanoplot", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("NanoPlot-report.html"), emit: report

   // conda 'bioconda::nanoplot=1.41.6'

    script:
    """
    NanoPlot --fastq ${reads} --maxlength 40000 --plots kde --legacy kde --no_static --outdir nanoplot_${sample_id}
    cp nanoplot_${sample_id}/* .
    """
}

// ==========================
// 3. MultiQC
// ==========================
process multiqc {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path(fastqc_reports)

    output:
    path("multiqc_report.html"), emit: report

  //  conda 'bioconda::multiqc=1.14'

    script:
    """
    mkdir -p fastqc_reports
    cp ${fastqc_reports}/* fastqc_reports/ 2>/dev/null || :
    multiqc fastqc_reports/ -n multiqc_report.html
    """
}

// ==========================
// 4. Flye assembly
// ==========================
process flye_assembly {
    tag "$sample_id"
    publishDir "${params.outdir}/assemblies", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("assembly.fasta"), emit: assembly

    conda '/home/user/anaconda3/envs/flye_env'

    script:
    """
    flye --nano-raw ${reads} --genome-size 2.2m --plasmids -o flye_${sample_id} --threads ${task.cpus}
    cp flye_${sample_id}/assembly.fasta .
    """
}

// ==========================
// 5. QUAST
// ==========================
process quast {
    tag "$sample_id"
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("quast_results/report.html"), emit: report

   // conda 'bioconda::quast=5.2.0'

    script:
    """
    quast.py -o quast_results ${assembly}
    """
}

// ==========================
// 6. Prokka
// ==========================
process prokka {
    tag "$sample_id"
    publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("prokka_output/*"), emit: annotation

  //  conda 'bioconda::prokka=1.14.6'

    script:
    """
    prokka --outdir prokka_output --prefix ${sample_id} --cpus ${task.cpus} --force ${assembly}
    """
}

// ==========================
// 7. Abricate
// ==========================
process abricate {
    tag "$sample_id"
    publishDir "${params.outdir}/abricate", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    path("*.tab"), emit: results

  //  conda 'bioconda::abricate=1.0.1'

    script:
    """
    abricate --db ncbi ${assembly} > resistance_genes.tab
    abricate --db plasmidfinder ${assembly} > plasmids.tab
    """
}

// ==========================
// Workflow
// ==========================
workflow {

    // 1. FastQC
    fastqc_out = reads_ch | fastqc

    // 2. NanoPlot
    nanoplot_out = reads_ch | nanoplot

    // 3. MultiQC after FastQC
    multiqc(fastqc_out.reports.collect())

    // 4. Flye assembly after MultiQC finishes
    flye_out = reads_ch | flye_assembly

    // 5. QUAST
    quast(flye_out.assembly)

    // 6. Prokka
    prokka(flye_out.assembly)

    // 7. Abricate
    abricate(flye_out.assembly)
}

// ==========================
// Completion message
// ==========================
workflow.onComplete {
    println """
    Pipeline completed successfully!

    Results directory: ${params.outdir}
    Assemblies: ${params.outdir}/assemblies
    Annotations: ${params.outdir}/annotations
    """
}

