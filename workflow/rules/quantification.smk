rule prepare_stu_nuclear_gtf:
    input:
        gff3=STU_GFF3
    output:
        gtf=STU_NUCLEAR_GTF
    params:
        script=SCRIPT_PATH + "/prepare_stu_gtf.py"
    conda:
        "../envs/quantification.yaml"
    log:
        WORK_DIR + "/00.reference/stu_nuclear_gtf.log"
    shell:
        """
        python {params.script:q} \
            --input {input.gff3:q} \
            --output {output.gtf:q} \
            > {log:q} 2>&1
        """


rule stu_nuclear_bam:
    input:
        bam=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.bam"
        ),
        bai=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/04.hisat2/{sample}.filtered.unique.bam.bai"
        )
    output:
        bam=(
            STU_NUCLEAR_BAM_DIR
            + "/{sample}/{sample}.stu_nuclear.bam"
        ),
        bai=(
            STU_NUCLEAR_BAM_DIR
            + "/{sample}/{sample}.stu_nuclear.bam.bai"
        ),
        contig_counts=(
            WORK_DIR
            + "/03.expression/{sample}/{sample}.stu_nuclear.contig_counts.tsv"
        )
    params:
        stu_nuclear_regex=config["reference"]["stu_nuclear_regex"],
        stu_organelle_regex=config["reference"]["stu_organelle_regex"]
    conda:
        "../envs/quantification.yaml"
    threads:
        config["quantification"]["samtools_threads"]
    resources:
        mem_mb=config["quantification"]["samtools_memory_mb"]
    log:
        WORK_DIR + "/03.expression/{sample}/{sample}.stu_nuclear_bam.log"
    shell:
        r"""
        mapfile -t nuclear_contigs < <(
            samtools idxstats {input.bam:q} \
            | awk \
                -v nuclear={params.stu_nuclear_regex:q} \
                -v organelle={params.stu_organelle_regex:q} \
                '$1 ~ nuclear && $1 !~ organelle {{print $1}}'
        )
        if (( ${{#nuclear_contigs[@]}} == 0 )); then
            echo "No Stu nuclear contigs found in BAM header" > {log:q}
            exit 1
        fi

        samtools view \
            -@ {threads} \
            -b \
            -o {output.bam:q} \
            {input.bam:q} \
            "${{nuclear_contigs[@]}}" \
            2> {log:q}
        samtools index -@ {threads} -o {output.bai:q} {output.bam:q}

        samtools idxstats {output.bam:q} \
        | awk 'BEGIN {{OFS="\t"; print "contig", "mapped_reads"}}
            $3 > 0 {{print $1, $3}}' \
        > {output.contig_counts:q}

        if awk \
            -v nuclear={params.stu_nuclear_regex:q} \
            -v organelle={params.stu_organelle_regex:q} \
            'NR > 1 && ($1 !~ nuclear || $1 ~ organelle) {{exit 1}}' \
            {output.contig_counts:q}; then
            :
        else
            echo "Non-nuclear alignments remain in {output.bam}" >> {log:q}
            exit 1
        fi
        """


rule stringtie_tpm:
    input:
        bam=(
            STU_NUCLEAR_BAM_DIR
            + "/{sample}/{sample}.stu_nuclear.bam"
        ),
        bai=(
            STU_NUCLEAR_BAM_DIR
            + "/{sample}/{sample}.stu_nuclear.bam.bai"
        ),
        gtf=STU_NUCLEAR_GTF,
        strand=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.strandness.tsv"
        )
    output:
        gtf=WORK_DIR + "/03.expression/{sample}/{sample}.stringtie.gtf",
        abundance=(
            WORK_DIR
            + "/03.expression/{sample}/{sample}.stringtie.gene_abund.tsv"
        )
    conda:
        "../envs/quantification.yaml"
    threads:
        config["quantification"]["stringtie_threads"]
    resources:
        mem_mb=config["quantification"]["stringtie_memory_mb"]
    log:
        WORK_DIR + "/03.expression/{sample}/{sample}.stringtie.log"
    shell:
        r"""
        strandness=$(awk 'NR == 2 {{print $6}}' {input.strand:q})
        strand_option=()
        if [[ "$strandness" == "FR" ]]; then
            strand_option=(--fr)
        elif [[ "$strandness" == "RF" ]]; then
            strand_option=(--rf)
        elif [[ "$strandness" != "unstranded" ]]; then
            echo "Invalid strandness value: $strandness" > {log:q}
            exit 1
        fi

        stringtie \
            {input.bam:q} \
            -p {threads} \
            -G {input.gtf:q} \
            -e \
            -A {output.abundance:q} \
            -o {output.gtf:q} \
            "${{strand_option[@]}}" \
            > {log:q} 2>&1
        """


rule featurecounts_readcount:
    input:
        bam=(
            STU_NUCLEAR_BAM_DIR
            + "/{sample}/{sample}.stu_nuclear.bam"
        ),
        bai=(
            STU_NUCLEAR_BAM_DIR
            + "/{sample}/{sample}.stu_nuclear.bam.bai"
        ),
        gtf=STU_NUCLEAR_GTF,
        strand=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/03.strand_analysis/{sample}.strandness.tsv"
        )
    output:
        counts=(
            WORK_DIR
            + "/03.expression/{sample}/{sample}.featurecounts.tsv"
        ),
        summary=(
            WORK_DIR
            + "/03.expression/{sample}/{sample}.featurecounts.tsv.summary"
        )
    conda:
        "../envs/quantification.yaml"
    threads:
        config["quantification"]["featurecounts_threads"]
    resources:
        mem_mb=config["quantification"]["featurecounts_memory_mb"]
    log:
        WORK_DIR + "/03.expression/{sample}/{sample}.featurecounts.log"
    shell:
        r"""
        strandness=$(awk 'NR == 2 {{print $6}}' {input.strand:q})
        if [[ "$strandness" == "FR" ]]; then
            strand_code=1
        elif [[ "$strandness" == "RF" ]]; then
            strand_code=2
        elif [[ "$strandness" == "unstranded" ]]; then
            strand_code=0
        else
            echo "Invalid strandness value: $strandness" > {log:q}
            exit 1
        fi

        featureCounts \
            -T {threads} \
            -p \
            --countReadPairs \
            -B \
            -C \
            -s "$strand_code" \
            -t exon \
            -g gene_id \
            -a {input.gtf:q} \
            -o {output.counts:q} \
            {input.bam:q} \
            > {log:q} 2>&1
        """


rule merge_quantification:
    input:
        tpm=expand(
            WORK_DIR
            + "/03.expression/{sample}/{sample}.stringtie.gene_abund.tsv",
            sample=SAMPLES
        ),
        counts=expand(
            WORK_DIR
            + "/03.expression/{sample}/{sample}.featurecounts.tsv",
            sample=SAMPLES
        )
    output:
        tpm=WORK_DIR + "/03.expression/matrices/gene_tpm.tsv",
        counts=WORK_DIR + "/03.expression/matrices/gene_readcount.tsv"
    params:
        script=SCRIPT_PATH + "/merge_quantification.py",
        samples=SAMPLES
    conda:
        "../envs/quantification.yaml"
    shell:
        """
        python {params.script:q} \
            --samples {params.samples:q} \
            --tpm-files {input.tpm:q} \
            --count-files {input.counts:q} \
            --output-tpm {output.tpm:q} \
            --output-count {output.counts:q}
        """
