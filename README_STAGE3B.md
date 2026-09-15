# Stage 3B 无标签 SAG 聚合（拖拉机版）

`sag-stage3b-tractor` 是独立 C++17 编排器，不修改现有
`dna2bit-sag-pipeline` 的 1–3A 逻辑。它负责输入闭合、外部程序调度、结果解析、
证据表、FASTA 合并和收据。SAG–SAG ANI/AF 计算路径接仓库内独立 C++17
`gtdb-ani-af` 的 exact triangle 接口；BLAST+、官方 `python-igraph` /
`leidenalg` 及已有 `cpp-subass-full` 仍是科学依赖。当前 3B **不再搜索
GTDB reference**，也不再使用 `--gtdb-taxonomy`、`--ani-db` 或
“double-negative”筛选。

**该接线是实验路径，不是生产默认。** `gtdb-ani-af` 不复制、链接或调用 skani，
但也不冒充 skani 数值模型。54 个历史 top-pair 的首轮 95/50 决策虽为 54/54
一致，ANI MAE 仍为 4.47，低相似区 AF 偏差大；G0045 因当前 ungapped 固定
diagonal 低估 ANI 4.89。运行必须显式给
`--allow-experimental-ani-engine`，真实生产替换尚未获准。

权威口径来自服务器：

- `/home/data/fyc/lake/25_recluster_1530_newmethod/global_signed.py`
- `/home/data/fyc/lake/25_recluster_1530_newmethod/marker_blastn_mp.py`
- `/home/data/fyc/lake/25_recluster_1530_newmethod/global_signed/report/pipeline_report_unpeeled.html`

网页入口为
`http://biotrainee.cn:13053/user/fyc/lab/tree/lake/25_recluster_1530_newmethod/global_signed/report/pipeline_report_unpeeled.html`。
软件、规范和测试统一使用该 `global_signed` 权威路径。

## 当前科学流程

1. 输入只能是老师原版 `dna2bit` 没有接受 species 标签的 SAG，并且已经通过
   更早的 assembly 总长度门 `total_bp >=1000`。3B 会严格核对 manifest 中的
   negative reason 和 FASTA 证据，不接受正标签 SAG 混入。
2. 在建图前，对全部 Dna2bit-negative SAG 应用 3B 质量门：
   `max_contig >=1000` 且 CheckM2 1.0.1 `contamination <5.0`。边界值
   `max_contig=1000` 通过，`contamination=5.0` 不通过。通过者就是完整图节点
   集合；这里不再搜索 GTDB reference，也不再减去 reference-positive SAG。
3. 对质量合格节点执行全部 `N*(N-1)/2` 个 SAG–SAG pair 的 exact triangle：

   ```bash
   gtdb-ani-af triangle --list cellbit_negative_quality_pass.list \
     --out ani_triangle_sparse.tsv --stats-out triangle_stats.json --threads N \
     --profile medium --edge-mode sparse --min-af 15 --triangle-mode exact
   ```

   `--triangle-mode exact` 表示所有节点 pair 都进入计算，不先靠 sketch 候选召回
   删除 pair；`--profile medium` 冻结本引擎 mapper 参数；`--edge-mode sparse
   --min-af 15` 只控制结果表保留 `max(AF_ref,AF_query)>=15` 的行。这些名字表达
   本引擎自己的透明语义，不宣称与旧 skani `-E --medium --min-af 15` 数值等价。
   后续仍仅保留 `ANI>=95` 且 `abs(GC差)<=2` 的吸引边：
   `w_pos=((ANI-95)/5)*(max(AF_ref,AF_query)/100)+0.01`。这里 AF 是权重，
   不是 pairwise 硬门。
4. 在 marker 比对前，吸引图先忠实执行一次历史 positive-only Leiden：
   `leidenalg.find_partition(RBConfigurationVertexPartition,
   resolution_parameter=1.0, seed=20260811, n_iterations=-1)`。只有其
   size `>=10` 社区成员进入 marker 阶段；最终 signed Leiden 仍覆盖全部质量
   合格节点，未做 marker 的节点只是没有负边。
5. bac120 核酸 marker map 中，每个目标 SAG/marker 只保留最长序列，且长度至少
   150 nt。每 marker 独立 `makeblastdb` / `blastn`，单 job 固定 2 线程：

   ```bash
   blastn -evalue 1e-5 -max_target_seqs 10000 \
     -outfmt '6 qseqid sseqid pident length'
   ```

   每个 SAG pair/marker 只取 alignment length 最长的一条；跨 marker 聚合。
   `n_markers>=3 && mean_pident<97` 才生成排斥边，
   `w_neg=(97-mean_pident)/97`。

   positive-only 预聚类没有 size>=10 目标，或 marker map 中没有长度>=150 且
   覆盖至少两个目标 SAG 的 marker，都是合法的“零负证据”情形：程序写出只有
   header 的 marker pair/negative edge 表并继续最终 signed Leiden，不会报错。
6. `python/stage3b_signed_leiden.py` 是参数冻结的 Python 后端，原因是 C++
   编排层没有依赖 Leiden C API。它直接调用官方 `igraph` / `leidenalg`：正层
   `RBConfigurationVertexPartition(resolution=r)`，负层
   `CPMVertexPartition(resolution=0)`，multiplex layer weights `[+1,-1]`，
   seed `20260811`，`n_iterations=-1`。扫描顺序为：
   `(1,1),(1,3),(1,10),(2,3),(2,10),(2,30)`。
   后端可显式指定 `--workers N`，将六组互不依赖的参数放入最多六个独立进程；
   每个进程仍各自设置同一冻结 seed，并完整运行 `n_iterations=-1`。评分、同分
   决胜、报告和 membership 始终由父进程按上述冻结顺序生成，因此单进程与
   多进程输出逐字节一致。后端单独调用时默认 `--workers 1`，保持历史串行行为；
   本速度优化版编排器传入 `--workers min(6, --threads)`，以少量额外内存并行
   六个参数点，不改变参数顺序、评分或同分决胜。
7. 对每次扫描，仅 size `>=10` 的簇参与汇总。簇 purity 是簇内所有
   `n_markers>=3` 的 marker pair 中 `mean_pident>=97` 的比例；无可评估 pair
   的簇不算 pure。按
   `(size>=10 且 purity>=0.9 的簇数, 这些 pure 簇的 SAG 总数)` 词典序最大
   选参，同分保持上述扫描顺序。**最终 membership 忠实输出选中参数下全部
   size>=10 簇**，并非只输出 pure 簇。
8. 每个输出簇合并所有成员 FASTA，header 精确改为
   `SAG_ID__原header`（双下划线），再调用已有 `cpp-subass-full`。簇表同时记录
   `n_sags` 与 `n_found`；本实现对缺失/空 FASTA fail-closed，因此正常闭包时二者
   必须相等。这是相对旧脚本“缺失/空 FASTA 跳过但记录”（例如历史 G0001
   `474/473`）的安全升级；默认不静默丢样本。历史底层等价于
   `flye --subassemblies INPUT --out-dir OUT --threads 16`，
   此处通过已验证的 `cpp-subass-full` 封装调用。

旧 GTDB-reference 双阴性流程曾报告选择 `r=2, lambda=30`、58 个 size>=10 簇。
这些数字只属于历史审计，不是当前“全部质量合格 Dna2bit-negative”入口的预期值。

## 必需的上游输入与未自动化边界

本版**不会**自动运行 CheckM2 1.0.1 或 GTDB-Tk identify，不能把 marker map
当作凭空存在。调用方必须显式提供：

- `--manifest`：现有 `03B_unclassified_pending.tsv`，列为
  `sag_id, assembly_fasta, reason`；`reason` 必须精确为
  `dna2bit_rejected_or_no_hit`、`dna2bit_negative`、`dna2bit_no_hit` 或
  `dna2bit_rejected` 之一。包含式模糊字符串不接受，因此
  `dna2bit_not_negative`、`positive_not_rejected` 会 fail-closed。
- `--quality-manifest`：至少含 `sag_id, max_contig, gc_pct,
  checkm2_contamination`，可含 `total_bp`。程序会重算 FASTA 总 bp、最大 contig
  和 GC，并拒绝不一致数据；GC 忠实使用 `(G+C)/全部序列字符`（N/ambiguity
  仍计入分母）并按旧表 4 位小数校验。CheckM2 contamination 是显式上游证据。
- `--marker-map`：至少含 `SAG_id, marker_id, sequence`。它应由 GTDB-Tk
  2.7.2 对同一批 SAG 执行下列 identify 命令后整理得到：

  ```bash
  gtdbtk identify --genome_dir SAG_DIR --out_dir GTDBTK_OUT -x fna \
    --cpus N --force --write_single_copy_genes
  ```

  历史整理脚本读取 `tigrfam`/`pfam` top-hit 表，按 bitscore 为每个
  SAG/marker 取 top hit，再从 `<SAG>_protein.fna` 提取对应核酸序列，形成至少
  `SAG_id, marker_id, sequence` 三列（也可保留 `gene_id, nt_len` 审计列）。这里
  GTDB-Tk 只用于取得 bac120 marker 序列，不参与 reference taxonomy 过滤。
- `--ani-engine`：自研可执行文件，默认 `gtdb-ani-af`；必须另给
  `--allow-experimental-ani-engine` 才允许执行。

因此，CheckM2 和 GTDB-Tk 的运行、版本收据以及从 single-copy genes 生成
`bac120_marker_nt_map.tsv` 仍是待接入缺口；缺少这些输入时程序会停止，不会降级。

## 编译、预演与运行

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

python3 -m pip install -r requirements-stage3b.txt

./build/sag-stage3b-tractor \
  --manifest output/03B_unclassified_pending.tsv \
  --quality-manifest stage3b_quality.tsv \
  --marker-map bac120_marker_nt_map.tsv \
  --ani-engine "$PWD/build/gtdb-ani-af" \
  --allow-experimental-ani-engine \
  --leiden-backend "$PWD/python/stage3b_signed_leiden.py" \
  --subass /path/to/cpp-subass-full \
  --out output/03B_unlabelled --threads 16 --dry-run
```

去掉 `--dry-run` 才会运行。中断后只能加 `--resume`。外部阶段使用
`sag-stage3b-receipt-v2`：恢复时重新解析 `status=PASS`，并校验 stage/command
SHA-256、显式输入集合指纹、config 指纹及每个输出的 size+SHA-256 汇总指纹。
内部确定性 TSV/FASTA 会从当前输入重算并逐字节/按 SHA-256 比对。旧收据、只凭
文件存在、无收据输出、混入另一轮的文件或被修改的输出均不能复用。

`--dry-run` 会验证所有输入数据与门槛，打印 exact triangle 计划，但不创建输出
或执行命令；后续 BLAST/Leiden/子装配依赖实际 triangle 结果，所以不会伪造
完整计划。

## 输出

- `quality_audit.tsv`：总 bp、max contig、GC、CheckM2 双门审计；
- `01_cellbit_negative_quality_pass/`：完整图节点表和 exact triangle 输入清单；
- `02_pairwise_ani_af/triangle_stats.json`：exact pair 数及输入/输出 SHA 闭包；
  `positive_edges.tsv` 是 GC 护栏后的加权吸引边；
- `03_positive_precluster/`：positive-only Leiden 的 marker 目标成员及收据；
- `04_marker_blastn/marker_pairwise_blastn.tsv`、`negative_edges.tsv`：marker 证据；
- `05_signed_leiden/chosen_membership.tsv`、`signed_report.json`：全部 size>=10
  簇与六组扫描审计；`unaggregated_sags.tsv` 显式列出孤立节点和小簇成员；
- `06_subassemble/<cluster>/`：双下划线 header 的合并 FASTA 和子装配结果；
- 每个外部阶段旁的 `*.PASS.json` 与顶层 `COMPLETE.json`。

顶层 `COMPLETE.json` 使用 `sag-stage3b-complete-v2`。它明确记录
`reference_search=disabled`、`reference_positive_exclusion=false`、
`triangle_mode=exact_all_pairs`、请求/实际评估的 pair 数，以及
`assigned_sags + unaggregated_sags == quality_pass_graph_nodes` 的入口—出口闭包。
`triangle_rows_emitted` 可以因 sparse `--min-af 15` 小于全 pair 数，但
`pairwise_pairs_evaluated` 必须严格等于 `N*(N-1)/2`。

## 微型测试

测试不启动 GTDB 全库、BLAST 或 Flye：

```bash
./build/sag-stage3b-tractor --self-test
tests/gtdb_ani_af_synthetic.sh ./build/gtdb-ani-af
python3 tests/run_stage3b_streaming_equivalence_smoke_20260904.py \
  --source src/stage3b.cpp
python3 python/stage3b_signed_leiden.py --self-test
python3 tests/test_stage3b_signed_leiden_parallel_equivalence_20260904.py
```

它们检查关键公式/accession 归一化、gzip/pair/index/search/triangle 小型接口，以及
purity、size>=10 全量输出与参数表。
