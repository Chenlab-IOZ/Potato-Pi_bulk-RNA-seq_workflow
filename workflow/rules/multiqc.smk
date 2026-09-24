rule multiqc:
    input:
        fastp=expand(
            WORK_DIR + "/01.pre_analysis/{sample}/01.cleandata/{sample}.fastp.json",
            sample=SAMPLES
        ),
        sortmerna=expand(
            WORK_DIR + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.sortmerna.log",
            sample=SAMPLES
        ),
        hisat2=expand(
            WORK_DIR + "/01.pre_analysis/{sample}/04.hisat2/{sample}.hisat2_summary.txt",
            sample=SAMPLES
        ),
        samtools=expand(
            WORK_DIR + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.flagstat.txt",
            sample=SAMPLES
        )
    output:
        html=WORK_DIR + f"/02.pre_assesment/01.multiQC/{MULTIQC_FILENAME}.html"
    params:
        outdir=WORK_DIR + "/02.pre_assesment/01.multiQC",
        filename=MULTIQC_FILENAME
    conda:
        "../envs/multiqc.yaml"
    threads:
        config["multiqc"]["threads"]
    resources:
        mem_mb=config["multiqc"]["memory_mb"]
    log:
        WORK_DIR + "/02.pre_assesment/01.multiQC/multiqc.log"
    shell:
        """
        multiqc {input} \
            --outdir {params.outdir} \
            --filename {params.filename} \
            > {log} 2>&1
        """
