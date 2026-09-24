rule sortmerna_index:
    input:
        reference=RRNA_REFERENCE
    output:
        index=directory(SORTMERNA_INDEX_DIR)
    params:
        workdir=WORK_DIR + "/00.reference/sortmerna_index"
    conda:
        "../envs/sortmerna.yaml"
    threads:
        config["sortmerna"]["index_threads"]
    resources:
        mem_mb=config["sortmerna"]["index_memory_mb"]
    log:
        WORK_DIR + "/00.reference/sortmerna_index/index.log"
    shell:
        """
        sortmerna \
            --ref {input.reference:q} \
            --workdir {params.workdir:q} \
            --threads {threads} \
            --index 1 \
            --task 5 \
            > {log:q} 2>&1
        """


rule sortmerna_filter:
    input:
        reads1=WORK_DIR + "/01.pre_analysis/{sample}/01.cleandata/{sample}_R1.clean.fastq.gz",
        reads2=WORK_DIR + "/01.pre_analysis/{sample}/01.cleandata/{sample}_R2.clean.fastq.gz",
        reference=RRNA_REFERENCE,
        index=rules.sortmerna_index.output.index
    output:
        non_rrna_r1=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA_fwd.fastq.gz"
        ),
        non_rrna_r2=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA_rev.fastq.gz"
        ),
        rrna_r1=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.rRNA_fwd.fastq.gz"
        ),
        rrna_r2=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.rRNA_rev.fastq.gz"
        ),
        stats=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.rRNA_stats.tsv"
        ),
        report=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.sortmerna.log"
        ),
        workdir=temp(directory(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/.work"
        ))
    params:
        aligned_prefix=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.rRNA"
        ),
        other_prefix=(
            WORK_DIR
            + "/01.pre_analysis/{sample}/02.sortmerna/{sample}.non_rRNA"
        )
    conda:
        "../envs/sortmerna.yaml"
    threads:
        config["sortmerna"]["threads"]
    resources:
        mem_mb=config["sortmerna"]["memory_mb"]
    log:
        WORK_DIR + "/01.pre_analysis/{sample}/logs/sortmerna.log"
    shell:
        """
        # SortMeRNA 7.0.0 indexed-gzip aborts on the full-size project FASTQs.
        # Keep the validated legacy feed until that upstream bug is fixed.
        sortmerna \
            --ref {input.reference:q} \
            --idx-dir {input.index:q} \
            --reads {input.reads1:q} \
            --reads {input.reads2:q} \
            --workdir {output.workdir:q} \
            --aligned {params.aligned_prefix:q} \
            --other {params.other_prefix:q} \
            --fastx \
            --paired_in \
            --out2 \
            --zip-out 1 \
            --readfeed 1 \
            --index 0 \
            --threads {threads} \
            > {log:q} 2>&1

        mv {params.other_prefix:q}_fwd.fq.gz {output.non_rrna_r1:q}
        mv {params.other_prefix:q}_rev.fq.gz {output.non_rrna_r2:q}
        mv {params.aligned_prefix:q}_fwd.fq.gz {output.rrna_r1:q}
        mv {params.aligned_prefix:q}_rev.fq.gz {output.rrna_r2:q}
        cp {log:q} {output.report:q}

        total_pairs=$(gzip -cd {input.reads1:q} | awk '
            END {{
                if (NR % 4 != 0) exit 1
                print NR / 4
            }}
        ')
        non_rrna_pairs=$(gzip -cd {output.non_rrna_r1:q} | awk '
            END {{
                if (NR % 4 != 0) exit 1
                print NR / 4
            }}
        ')
        rrna_pairs=$((total_pairs - non_rrna_pairs))

        awk \
            -v sample="{wildcards.sample}" \
            -v total="$total_pairs" \
            -v non_rrna="$non_rrna_pairs" \
            -v rrna="$rrna_pairs" \
            'BEGIN {{
                OFS="\t"
                fraction=(total > 0 ? rrna / total : 0)
                print "sample", "total_pairs", "non_rrna_pairs", "rrna_pairs", "rrna_fraction"
                print sample, total, non_rrna, rrna, fraction
            }}' \
            > {output.stats:q}
        """
