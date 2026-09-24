rule prepare_strand_reference_bed:
    input:
        stu_bed=STU_GENE_BED,
        pi_gtf=PI_GTF
    output:
        bed=COMBINED_GENE_BED
    shell:
        r"""
        (
            awk 'BEGIN {{OFS="\t"}}
                {{
                    chrom=$1
                    if (chrom !~ /^Stu_/) chrom="Stu_" chrom
                    start=$2
                    end=$3
                    size=end-start
                    print chrom,start,end,$4,0,$6,start,end,0,1,size,0
                }}' {input.stu_bed:q}

            awk -F '\t' 'BEGIN {{OFS="\t"}}
                $3 == "gene" {{
                    chrom=$1
                    if (chrom !~ /^Pinf_/) chrom="Pinf_" chrom
                    start=$4-1
                    end=$5
                    name=$9
                    sub(/.*gene_id "/, "", name)
                    sub(/".*/, "", name)
                    size=end-start
                    print chrom,start,end,name,0,$7,start,end,0,1,size,0
                }}' {input.pi_gtf:q}
        ) | LC_ALL=C sort -k1,1 -k2,2n > {output.bed:q}
        """


rule strand_subsample:
    input:
        reads1=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA_fwd.fastq.gz"
        ),
        reads2=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA_rev.fastq.gz"
        )
    output:
        reads1=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.subsample_R1.fastq.gz"
        ),
        reads2=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.subsample_R2.fastq.gz"
        )
    params:
        pairs=config["strand_analysis"]["sample_pairs"],
        seed=config["strand_analysis"]["random_seed"]
    conda:
        "../envs/alignment.yaml"
    shell:
        """
        seqtk sample -s{params.seed} {input.reads1:q} {params.pairs} \
            | gzip -c > {output.reads1:q}
        seqtk sample -s{params.seed} {input.reads2:q} {params.pairs} \
            | gzip -c > {output.reads2:q}
        """


rule strand_prealign:
    input:
        reads1=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.subsample_R1.fastq.gz"
        ),
        reads2=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.subsample_R2.fastq.gz"
        ),
        index=HISAT2_INDEX_FILES
    output:
        bam=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.prealign.bam"
        ),
        bai=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.prealign.bam.bai"
        )
    params:
        index_prefix=HISAT2_INDEX_PREFIX,
        max_alignments=config["strand_analysis"]["max_alignments"],
        max_intron_length=config["hisat2"]["max_intron_length"],
        align_threads=config["strand_analysis"]["hisat2_threads"],
        view_threads=config["strand_analysis"]["samtools_view_threads"],
        sort_threads=config["strand_analysis"]["samtools_sort_threads"]
    conda:
        "../envs/alignment.yaml"
    threads:
        config["strand_analysis"]["threads"]
    resources:
        mem_mb=config["strand_analysis"]["memory_mb"]
    log:
        WORK_DIR + "/01.pre_analysis/{sample}/logs/strand_prealign.log"
    shell:
        """
        hisat2 \
            -p {params.align_threads} \
            -k {params.max_alignments} \
            --max-intronlen {params.max_intron_length} \
            --no-mixed \
            --no-discordant \
            -x {params.index_prefix:q} \
            -1 {input.reads1:q} \
            -2 {input.reads2:q} \
            2> {log:q} \
        | samtools view -@ {params.view_threads} -b -f 2 -F 2304 - \
        | samtools sort -@ {params.sort_threads} -o {output.bam:q} -

        samtools index -@ {params.sort_threads} -o {output.bai:q} {output.bam:q}
        """


rule infer_strandedness:
    input:
        bam=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.prealign.bam"
        ),
        bai=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.prealign.bam.bai"
        ),
        bed=COMBINED_GENE_BED
    output:
        raw=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.infer_experiment.txt"
        ),
        tsv=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.strandness.tsv"
        )
    params:
        parser=SCRIPT_PATH + "/parse_strandness.py",
        sample_size=config["strand_analysis"]["sample_pairs"],
        stranded_threshold=config["strand_analysis"]["stranded_threshold"],
        unstranded_threshold=config["strand_analysis"]["unstranded_threshold"]
    conda:
        "../envs/rseqc.yaml"
    log:
        WORK_DIR + "/01.pre_analysis/{sample}/logs/infer_strandedness.log"
    shell:
        """
        infer_experiment.py \
            -i {input.bam:q} \
            -r {input.bed:q} \
            -s {params.sample_size} \
            > {output.raw:q} \
            2> {log:q}

        python {params.parser:q} \
            --input {output.raw:q} \
            --output {output.tsv:q} \
            --sample {wildcards.sample:q} \
            --stranded-threshold {params.stranded_threshold} \
            --unstranded-threshold {params.unstranded_threshold}
        """


rule merge_strandness:
    input:
        expand(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.strandness.tsv",
            sample=SAMPLES
        )
    output:
        WORK_DIR + "/02.pre_assesment/02.strandness/strandness_summary.tsv"
    shell:
        """
        awk 'FNR == 1 && NR != 1 {{next}} {{print}}' {input:q} > {output:q}
        """
