rule hisat2_index:
    input:
        genome=COMBINED_GENOME
    output:
        HISAT2_INDEX_FILES
    params:
        prefix=HISAT2_INDEX_PREFIX
    conda:
        "../envs/alignment.yaml"
    threads:
        config["hisat2"]["index_threads"]
    resources:
        mem_mb=config["hisat2"]["index_memory_mb"]
    log:
        WORK_DIR + "/00.reference/hisat2_index/index.log"
    shell:
        """
        hisat2-build \
            -p {threads} \
            {input.genome:q} \
            {params.prefix:q} \
            > {log:q} 2>&1
        """


rule hisat2_align_filter:
    input:
        reads1=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA_fwd.fastq.gz"
        ),
        reads2=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA_rev.fastq.gz"
        ),
        strand=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.strandness.tsv"
        ),
        index=HISAT2_INDEX_FILES
    output:
        bam=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.bam"
        ),
        bai=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.bam.bai"
        ),
        summary=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.hisat2_summary.txt"
        ),
        flagstat=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.flagstat.txt"
        )
    params:
        index_prefix=HISAT2_INDEX_PREFIX,
        max_alignments=config["hisat2"]["max_alignments"],
        max_intron_length=config["hisat2"]["max_intron_length"],
        minimum_mapq=config["hisat2"]["minimum_mapq"],
        align_threads=config["hisat2"]["align_threads"],
        view_threads=config["hisat2"]["samtools_view_threads"],
        sort_threads=config["hisat2"]["samtools_sort_threads"],
        post_threads=config["hisat2"]["samtools_post_threads"]
    conda:
        "../envs/alignment.yaml"
    threads:
        config["hisat2"]["threads"]
    resources:
        mem_mb=config["hisat2"]["memory_mb"]
    log:
        WORK_DIR + "/01.pre_analysis/{sample}/logs/hisat2.log"
    shell:
        r"""
        strandness=$(awk 'NR == 2 {{print $6}}' {input.strand:q})
        strand_option=()
        if [[ "$strandness" == "FR" || "$strandness" == "RF" ]]; then
            strand_option=(--rna-strandness "$strandness")
        elif [[ "$strandness" != "unstranded" ]]; then
            echo "Invalid strandness value: $strandness" >&2
            exit 1
        fi

        (
            hisat2 \
                -p {params.align_threads} \
                -k {params.max_alignments} \
                --max-intronlen {params.max_intron_length} \
                --no-mixed \
                --no-discordant \
                --new-summary \
                --summary-file {output.summary:q} \
                "${{strand_option[@]}}" \
                -x {params.index_prefix:q} \
                -1 {input.reads1:q} \
                -2 {input.reads2:q} \
            | awk 'BEGIN {{OFS="\t"}} /^@/ || /(^|\t)NH:i:1(\t|$)/' \
            | samtools view \
                -@ {params.view_threads} \
                -b \
                -f 3 \
                -F 3852 \
                -q {params.minimum_mapq} \
                - \
            | samtools sort \
                -@ {params.sort_threads} \
                -o {output.bam:q} \
                -
        ) 2> {log:q}

        samtools index -@ {params.post_threads} -o {output.bai:q} {output.bam:q}
        samtools flagstat -@ {params.post_threads} {output.bam:q} > {output.flagstat:q}
        """


rule mapping_composition:
    input:
        reads1=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA_fwd.fastq.gz"
        ),
        bam=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.bam"
        ),
        bai=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.bam.bai"
        )
    output:
        WORK_DIR
        + "/01.pre_analysis/{sample}/04.hisat2/{sample}.mapping_composition.tsv"
    params:
        stu_nuclear_regex=config["reference"]["stu_nuclear_regex"],
        stu_organelle_regex=config["reference"]["stu_organelle_regex"],
        pi_regex=config["reference"]["pi_regex"]
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        input_pairs=$(gzip -cd {input.reads1:q} | awk '
            END {{
                if (NR % 4 != 0) exit 1
                print NR / 4
            }}
        ')

        samtools view -f 64 {input.bam:q} \
        | awk \
            -v sample={wildcards.sample:q} \
            -v input_pairs="$input_pairs" \
            -v stu_nuclear_regex={params.stu_nuclear_regex:q} \
            -v stu_organelle_regex={params.stu_organelle_regex:q} \
            -v pi_regex={params.pi_regex:q} \
            'BEGIN {{OFS="\t"}}
            {{
                total++
                if ($3 ~ stu_organelle_regex) organelle++
                else if ($3 ~ stu_nuclear_regex) stu++
                else if ($3 ~ pi_regex) pi++
                else other++
            }}
            END {{
                mapped_fraction=(input_pairs > 0 ? total / input_pairs : 0)
                stu_fraction=(total > 0 ? stu / total : 0)
                organelle_fraction=(total > 0 ? organelle / total : 0)
                pi_fraction=(total > 0 ? pi / total : 0)
                other_fraction=(total > 0 ? other / total : 0)
                print "sample", "input_non_rrna_pairs", "filtered_mapped_pairs", \
                    "filtered_mapping_fraction", "Stu_nuclear_pairs", \
                    "Stu_nuclear_fraction", "MtPt_pairs", "MtPt_fraction", \
                    "Pi_pairs", "Pi_fraction", "other_pairs", "other_fraction"
                print sample, input_pairs, total, mapped_fraction, stu, \
                    stu_fraction, organelle, organelle_fraction, pi, \
                    pi_fraction, other, other_fraction
            }}' > {output:q}
        """


rule merge_mapping_composition:
    input:
        expand(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.mapping_composition.tsv",
            sample=SAMPLES
        )
    output:
        WORK_DIR + "/02.pre_assesment/03.mapping/mapping_composition.tsv"
    shell:
        """
        awk 'FNR == 1 && NR != 1 {{next}} {{print}}' {input:q} > {output:q}
        """
