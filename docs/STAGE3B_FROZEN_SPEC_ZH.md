# 3B：DNA2bit 无标签 SAG 的符号图聚类与共组装

> **历史冻结证据，已被当前生产规格取代。** 本文只保存老师旧
> `global_signed`/GTDB-reference 双阴性流程的可追溯规则，不是当前 Stage3B 的
> 运行说明。当前入口见 `README_STAGE3B.md`：不做 GTDB reference 过滤，全部
> 质量合格的 Dna2bit-negative SAG 进入 exact SAG–SAG 聚类。

本文件冻结老师现有 `global_signed` 流程的科学规则，防止封装时把网页摘要误当成完整实现。

权威来源：

- 报告：`/home/data/fyc/lake/25_recluster_1530_newmethod/global_signed/report/pipeline_report_unpeeled.html`
- 符号聚类：`/home/data/fyc/lake/25_recluster_1530_newmethod/global_signed.py`
- marker 比对：`/home/data/fyc/lake/25_recluster_1530_newmethod/marker_blastn_mp.py`
- 旧严格入口：`/home/data/fyc/lake/80_workflow_code/run_lake_cellbit_unknown_ppt_cluster_strict.sh`
- marker 提取：`/home/data/fyc/lake/80_workflow_code/prepare_lake_unknown_strict_cluster_inputs.py`

## 输入边界

3B 只接收满足以下条件的 SAG assembly：

1. 已通过本软件新增的 SAG 总碱基数硬门：`total_bp >= 1000`；
2. 老师原版 DNA2bit 没有接受 species 标签；
3. GTDB R232 skani best hit 也没有达到物种级门；
4. 旧 3B 质量门仍单独保留：CheckM2 `contamination < 5.0` 且 `max_contig >= 1000 bp`；
5. GC 可以从 assembly 计算。

第 1 条是本软件新增的前置硬门，不能冒充或替代第 4 条。

## 3B-1：skani 已知物种排除

运行：

```text
skani search -d <GTDB_R232_SKANI_DB> --ql <query.list> -t <threads> -o <search.tsv>
```

每个 SAG 在全部命中中按 `(ANI, AF_query, AF_ref)` 从大到小选唯一 best hit。参考 accession 必须能映射到非空 GTDB taxonomy。只有同时满足：

- `ANI >= 95.0`
- `AF_query >= 50.0`

才判为 skani species-positive，并从 3B 中排除。注意这里必须使用 `AF_query`，不能换成 `AF_ref`、两者平均值或 `max(AF_ref, AF_query)`。

DNA2bit-negative 且 skani-negative 的 SAG 才叫 double-negative/neither，并进入后续图聚类。

## 3B-2：质量门与 bac120 marker

旧 3B 使用 CheckM2 1.0.1 对单 SAG 做质量估计，保留：

- `contamination < 5.0`（严格小于 5，不包含 5.0）；
- assembly 的 `max_contig >= 1000 bp`。

对保留 SAG 运行：

```text
gtdbtk identify --genome_dir <sag_fasta_dir> --out_dir <gtdbtk_out> \
  -x fna --cpus <threads> --force --write_single_copy_genes
```

从 GTDB-Tk `tigrfam` 和 `pfam` top-hit 表中，按 bitscore 为每个 SAG 的每个 marker 选择 top hit，再从 `<SAG>_protein.fna` 取相应核酸序列，形成：

```text
SAG_id  marker_id  gene_id  nt_len  sequence
```

## 3B-3：吸引边

对 double-negative 质量合格 SAG 运行：

```text
skani triangle -l <assembly.list> -o <triangle.tsv> \
  -t <threads> -E --medium --min-af 15
```

候选对只有同时满足以下条件才形成吸引边：

- `ANI >= 95.0`
- `abs(GC_A - GC_B) <= 2.0` 个百分点

AF 不再设 50% 硬门，而取 `AF = max(AF_ref, AF_query)` 作为边权的一部分：

```text
w_pos = ((ANI - 95.0) / 5.0) * (AF / 100.0) + 0.01
```

同一对若出现多次，只保留 `w_pos` 最大的一条边。所有保留 SAG 都必须成为节点，包括没有任何边的孤立节点。

## 3B-4：marker 计算目标预聚类

原流程没有对全部 double-negative 节点运行昂贵的 marker all-vs-all。它先只用完整吸引图做一次普通 Leiden：

- 调用 `leidenalg.find_partition`
- partition 类型为 `RBConfigurationVertexPartition`
- `resolution_parameter = 1.0`
- RNG seed `20260811`
- `n_iterations = -1`

这里应忠实使用旧 `pipeline_global.py` 的 `find_partition` 接口；它没有像最终 signed Leiden 那样显式覆盖 `consider_comms`。

只有这次正图预聚类中 `size>=10` 社区的成员进入 bac120 marker 比对。历史数据因此从 7,408 个最终图节点缩到 3,963 个 marker 目标。这个预聚类只决定“给谁计算负证据”，不决定最终输出；最终符号 Leiden 仍须包含全部 double-negative 节点，其余节点只是没有排斥边。

## 3B-5：bac120 排斥边

预处理规则：

- marker 核酸序列长度必须 `>=150 nt`；
- 同一 SAG 同一 marker 若有多条序列，只保留最长一条；
- 一个 marker 至少出现在两个 SAG 中才运行 all-vs-all。

每个 marker 独立运行：

```text
makeblastdb -in <marker.fna> -dbtype nucl -out <marker_db>
blastn -query <marker.fna> -db <marker_db> \
  -outfmt "6 qseqid sseqid pident length" \
  -max_target_seqs 10000 -evalue 1e-5 \
  -num_threads 2 -out <marker.tsv>
```

旧运行把 marker 作业并发数设为 96、每作业固定 2 个 blastn 线程；封装允许按机器资源调小外层并发，但不能把单作业的 `-num_threads 2` 偷换成全局线程数。

对同一 SAG 对、同一 marker 的多条 alignment，只取 alignment length 最长的一条及其 pident。随后跨 marker 计算 `mean_pident`。仅当：

- `n_markers >= 3`
- `mean_pident < 97.0`

才建立排斥边：

```text
w_neg = (97.0 - mean_pident) / 97.0
```

## 3B-6：符号 Leiden

节点按 SAG ID 排序，构造两个同节点集的图层：

- 正层：`RBConfigurationVertexPartition`，权重 `w_pos`，分辨率 `r`；
- 负层：先把 `w_neg` 乘 `lambda`，再用 `CPMVertexPartition`，`resolution_parameter=0.0`；
- 多层权重固定为 `[+1, -1]`。

严禁把负层改回 RBConfiguration；其零模型项在负层符号下会产生错误方向。

优化器固定：

- RNG seed：`20260811`
- `consider_comms = ALL_COMMS`
- `n_iterations = -1`（运行至收敛）

只扫描六组参数：

```text
(r, lambda) =
(1.0, 1.0), (1.0, 3.0), (1.0, 10.0),
(2.0, 3.0), (2.0, 10.0), (2.0, 30.0)
```

每个簇的 marker purity 定义为：

```text
簇内 n_markers>=3 且 mean_pident>=97 的 SAG 对数量
---------------------------------------------------
簇内所有 n_markers>=3 的 SAG 对数量
```

参数选择分数为二元组：

```text
(size>=10 且 purity>=0.9 的簇数,
 这些纯簇包含的 SAG 总数)
```

按 Python 二元组的词典序取最大值；完全同分时保留扫描顺序中先出现的组合。历史湖泊数据由此选中 `r=2.0, lambda=30`，但新数据必须重新走同一六点扫描，不能把历史最优参数硬编码成普适最优。

最终 membership 写出该最优划分中的 **全部 `size>=10` 簇**。这里不能误改为“只写 purity>=0.9 的簇”；purity 门用于选参数，不是最终写簇过滤器。

## 3B-7：按簇共组装

把每个簇内所有 SAG assembly 的 contig 合并；每个 contig header 精确改为 `SAG_ID__原header`（两个下划线）以保留来源，序列不改。不存在或为空的 assembly 不得伪造成有效输入，且必须在 manifest 中同时记录簇的 `n_sags` 与实际读到的 `n_found`。随后调用本项目已有的 C++ Flye-subassemblies 全流程实现，科学模式等价于：

```text
flye --subassemblies <combined.fasta> --out-dir <cluster_out> --threads <threads>
```

只有 `size>=10` 的簇进入共组装。孤立节点、小于 10 个 SAG 的簇及未进入任何输出簇的 SAG 必须写入单独的未聚合清单，不能静默丢失。

## 明确不属于 3B 软件本体的步骤

报告第 8 节之后的 pathway-constrained 去冗余、冲突剥离、候选 CheckM2/GUNC 和 MIMAG 终审属于下游评估/精修，不并入当前 3B 聚合+subassemble 主体。尤其报告已明确指出 pathway 步骤会直接优化 CheckM2 指标且不等于生物学去嵌合。
