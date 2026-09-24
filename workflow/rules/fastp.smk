rule fastp:
    input:
        reads1=lambda wildcards: sample_fastq(wildcards, "r1"),
        reads2=lambda wildcards: sample_fastq(wildcards, "r2")
    output:
        creads1=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/01.cleandata/{sample}_R1.clean.fastq.gz"
        ),
        creads2=temp(
            WORK_DIR
            + "/01.pre_analysis/{sample}/01.cleandata/{sample}_R2.clean.fastq.gz"
        ),
        json=WORK_DIR + "/01.pre_analysis/{sample}/01.cleandata/{sample}.fastp.json",
        html=WORK_DIR + "/01.pre_analysis/{sample}/01.cleandata/{sample}.fastp.html"
    params:
        base_phred=SCRIPT_PATH + "/base_phred.sh",
        n_base_limit=config["fastp"]["n_base_limit"],
        average_qual=config["fastp"]["average_qual"],
        cut_front="--cut_front" if config["fastp"]["cut_front"] else "",
        cut_tail="--cut_tail" if config["fastp"]["cut_tail"] else "",
        cut_right="--cut_right" if config["fastp"]["cut_right"] else "",
        detect_adapter=(
            "--detect_adapter_for_pe"
            if config["fastp"]["detect_adapter_for_pe"]
            else ""
        )
    conda:
        "../envs/fastp.yaml"
    threads:
        config["fastp"]["threads"]
    resources:
        mem_mb=config["fastp"]["memory_mb"]
    log:
        WORK_DIR + "/01.pre_analysis/{sample}/logs/fastp.log"
    shell:
        """
        phred_value=$(bash {params.base_phred} {input.reads1})

        if [[ "$phred_value" == "64" ]]; then
            fastp \
                --thread {threads} \
                --phred64 \
                --in1 {input.reads1} \
                --in2 {input.reads2} \
                --out1 {output.creads1} \
                --out2 {output.creads2} \
                --json {output.json} \
                --html {output.html} \
                {params.cut_front} \
                {params.cut_tail} \
                {params.cut_right} \
                --n_base_limit {params.n_base_limit} \
                --average_qual {params.average_qual} \
                {params.detect_adapter} \
                --report_title "{wildcards.sample}" \
                > {log} 2>&1
        else
            fastp \
                --thread {threads} \
                --in1 {input.reads1} \
                --in2 {input.reads2} \
                --out1 {output.creads1} \
                --out2 {output.creads2} \
                --json {output.json} \
                --html {output.html} \
                {params.cut_front} \
                {params.cut_tail} \
                {params.cut_right} \
                --n_base_limit {params.n_base_limit} \
                --average_qual {params.average_qual} \
                {params.detect_adapter} \
                --report_title "{wildcards.sample}" \
                > {log} 2>&1
        fi
        """
