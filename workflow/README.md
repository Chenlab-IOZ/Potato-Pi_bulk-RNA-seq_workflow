# sample_relationship

在 `config/config.yaml` 中选择分析方法：

```yaml
sample_relationship:
  method: "both"  # deseq2 / edgeR / both
  min_total_count: 10
  top_variable_genes: 500
  correlation_method: "pearson"
  clustering_method: "average"
  memory_mb: 8192
```

当前默认 `both`；旧配置未指定 `method` 时使用 `edgeR`。无效选项会报错。
仅运行所选的样本关系分析（在项目根目录执行）：

```bash
snakemake -s 03_script/workflow/Snakefile sample_relationship_all \
  --cores 2 --use-conda --conda-prefix 03_script/.snakemake/conda
```

完整流程的 `all` 目标也包含所选方法。流程入口仍检查样本表、原始 FASTQ 和参考文件路径。
每种方法单独运行，结果写入配置的 `output_dir/04.sample_relationship/deseq2/`
或 `output_dir/04.sample_relationship/edgeR/`。切换方法不会删除已有结果。
原先直接位于 `04.sample_relationship/` 下的历史文件保留，但不再作为新流程的输出目标。

| 文件 | deseq2 | edgeR |
| --- | --- | --- |
| `normalized_expression.tsv` | 原始计数除以 DESeq2 size factor | TMM 有效文库大小标准化的 CPM（非 log） |
| `transformed_expression.tsv` | rlog（`blind=TRUE`） | log2-CPM（`prior.count=2`） |
| `normalization_factors.tsv` | 原始文库大小、size factor | 原始文库大小、TMM factor、有效文库大小 |
| `pca_scores.tsv` / `pca.pdf` | rlog 的 PCA 得分及图 | log2-CPM 的 PCA 得分及图 |
| `sample_correlation.tsv` / `sample_correlation_heatmap.pdf` | rlog 的相关性矩阵及热图 | log2-CPM 的相关性矩阵及热图 |
| `sample_clustering.pdf` | rlog 的欧氏距离层次聚类 | log2-CPM 的欧氏距离层次聚类 |
| `sample_relationship.log` | 过滤统计、方法及 R sessionInfo | 过滤统计、方法及 R sessionInfo |

两种方法均保留总计数严格大于 `min_total_count` 的基因。DESeq2 使用 median-ratio
size factors（`type="ratio"`）、`design=~1` 和盲 rlog 转换；不进行差异表达检验或批次校正。
如果没有跨所有样本均为正计数的基因，DESeq2 的 ratio 估计会报错，不自动更换算法。
PCA 对全部保留基因进行中心化，不缩放。相关性在各方法转换后分别按 MAD
选择最高的 `top_variable_genes` 个非零 MAD 基因，因此两个方法的相关性基因集可能不同。
热图使用 `1 - correlation` 距离；样本聚类使用全部保留基因的欧氏距离。
`normalized_expression.tsv` 在两种方法中的单位不同，不宜直接比较数值。

两种方法使用相同文件名便于下游读取；edgeR 的 `transformed_expression.tsv`
对应旧版 `logcpm_expression.tsv`。外部作图脚本若仍指向旧路径，需显式选择新方法目录。
