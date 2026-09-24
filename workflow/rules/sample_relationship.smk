# Older configurations without method retain the existing edgeR behavior.
RELATIONSHIP_METHOD = config["sample_relationship"].get("method", "edgeR")
if RELATIONSHIP_METHOD not in ("deseq2", "edgeR", "both"):
    raise WorkflowError("sample_relationship.method must be deseq2, edgeR, or both")
RELATIONSHIP_METHODS = (
    ["deseq2", "edgeR"] if RELATIONSHIP_METHOD == "both"
    else [RELATIONSHIP_METHOD]
)
RELATIONSHIP_FILES = [
    "normalized_expression.tsv", "transformed_expression.tsv",
    "normalization_factors.tsv", "sample_correlation.tsv", "pca_scores.tsv",
    "sample_correlation_heatmap.pdf", "pca.pdf", "sample_clustering.pdf",
]
SAMPLE_RELATIONSHIP_OUTPUTS = expand(
    WORK_DIR + "/04.sample_relationship/{method}/{filename}",
    method=RELATIONSHIP_METHODS, filename=RELATIONSHIP_FILES,
)


rule sample_relationship_all:
    input:
        SAMPLE_RELATIONSHIP_OUTPUTS


rule sample_relationship:
    input:
        SAMPLE_RELATIONSHIP_OUTPUTS


rule sample_relationship_method:
    input:
        counts=WORK_DIR + "/03.expression/matrices/gene_readcount.tsv",
        samples=str(SAMPLE_FILE),
        script=SCRIPT_PATH + "/sample_relationship.R"
    output:
        normalized=WORK_DIR + "/04.sample_relationship/{method}/normalized_expression.tsv",
        transformed=WORK_DIR + "/04.sample_relationship/{method}/transformed_expression.tsv",
        factors=WORK_DIR + "/04.sample_relationship/{method}/normalization_factors.tsv",
        correlation=WORK_DIR + "/04.sample_relationship/{method}/sample_correlation.tsv",
        pca_scores=WORK_DIR + "/04.sample_relationship/{method}/pca_scores.tsv",
        correlation_heatmap=(
            WORK_DIR
            + "/04.sample_relationship/{method}/sample_correlation_heatmap.pdf"
        ),
        pca=WORK_DIR + "/04.sample_relationship/{method}/pca.pdf",
        clustering=WORK_DIR + "/04.sample_relationship/{method}/sample_clustering.pdf"
    wildcard_constraints:
        method="|".join(RELATIONSHIP_METHODS)
    params:
        min_total_count=config["sample_relationship"]["min_total_count"],
        top_variable_genes=config["sample_relationship"]["top_variable_genes"],
        correlation_method=config["sample_relationship"]["correlation_method"],
        clustering_method=config["sample_relationship"]["clustering_method"]
    conda:
        "../envs/sample_relationship.yaml"
    resources:
        mem_mb=config["sample_relationship"]["memory_mb"]
    log:
        WORK_DIR + "/04.sample_relationship/{method}/sample_relationship.log"
    shell:
        """
        Rscript {input.script:q} \
            {input.counts:q} \
            {input.samples:q} \
            {output.transformed:q} \
            {output.correlation:q} \
            {output.pca_scores:q} \
            {output.correlation_heatmap:q} \
            {output.pca:q} \
            {output.clustering:q} \
            {params.min_total_count:q} \
            {params.top_variable_genes:q} \
            {params.correlation_method:q} \
            {params.clustering_method:q} \
            {wildcards.method:q} \
            {output.normalized:q} \
            {output.factors:q} \
            > {log:q} 2>&1
        """
