# `gtdb-ani-af`：独立 C++17 ANI/AF 估计器

`gtdb-ani-af` 是本仓库自行实现、算法完全可解释的 FASTA 相似度估计器。它没有
复制、链接或调用 skani，也没有训练模型、隐藏参数或 GTDB 私有算法。名字中的
GTDB 仅表示它能给 GTDB reference manifest 建索引；算法本身适用于普通核酸
FASTA。

**它不是 skani 的数值兼容实现，也不能把输出冒充 skani 结果。** 当前 Stage 3B
只接入它的 exact SAG–SAG `triangle` 路径，不使用 GTDB reference `index/search`
来过滤节点。数据库接口仍是独立实验能力；没有完成真实 SAG×GTDB 的召回率、
ANI/AF 偏差、阈值翻转率和资源对照前，不能宣称它可替代 skani reference search。

## 接口

单对比较：

```bash
gtdb-ani-af pair --query query.fna --reference reference.fna \
  --out pair.tsv
```

建立数据库候选索引。manifest 可为每行一个 FASTA 路径，或带 `Ref_file`、
`path`、`reference` 表头的第一列 TSV；相对路径以 manifest 所在目录解析。
taxonomy、species-radii 和 release metadata 在本实验引擎里**只作为版本来源证据**，
不参与 ANI/AF 计算或物种判定，但为了防止数据库和注释版本脱钩，建库时必须显式
提供并绑定 SHA-256：

```bash
gtdb-ani-af index --manifest gtdb_refs.tsv \
  --taxonomy bac120_taxonomy.tsv --radii species_radii.tsv \
  --release-metadata gtdb_release.txt --out-dir gtdb_ani_index \
  --k 15 --sketch-scale 64 --query-sketch-size 1000 \
  --max-posting 1000 --threads 8
```

批量 search：

```bash
gtdb-ani-af search --index gtdb_ani_index --queries queries.list \
  --out search.tsv --max-candidates 20 --top 5 --threads 8
```

若一个已经验证为 PASS 的索引只绑定错了 taxonomy release，不要重建 ANI/AF
核心索引。使用显式、write-once 的 rebind：

```bash
gtdb-ani-af rebind-taxonomy \
  --index gtdb_ani_index_v4_wrong_taxonomy \
  --taxonomy /home/data/shared/software/Dna2bit/GTDB232/genome_taxonomy_1.csv \
  --out-dir gtdb_ani_index_v4_gtdb232_taxonomy
```

命令会完整加载并验证源 PASS 索引，将新 taxonomy accession 集合与 `REFS.tsv`
归一化后的 accession 集合做严格全等比较，并拒绝重复、缺失或额外记录。新目录内的
`REFS.tsv`、`SKETCHES.bin`、`POSTINGS.bin`（以及 v4 的
`POSTINGS_LOOKUP.bin`）逐文件 SHA 必须与源收据完全相同；实现会尝试 hardlink，失败时
自动 copy，但从不相信文件操作本身，发布前总会重新计算 SHA。新目录还保存原始
`SOURCE_COMPLETE.json` 和新 taxonomy 的原字节副本 `TAXONOMY.tsv`，最终
`COMPLETE.json` 绑定源 index/receipt SHA、新 taxonomy 路径/SHA、accession 集合 SHA
及实际 materialization 模式。目标目录已存在时拒绝执行。

同一样本内部的 triangle 默认仍运行全部 `N(N-1)/2` 对；`exact` 是默认值，兼容
旧命令和旧输出顺序：

```bash
gtdb-ani-af triangle --list assemblies.list --out triangle.tsv \
  --stats-out triangle.stats.json \
  --threads 96 --profile medium --edge-mode sparse --min-af 15 \
  --triangle-mode exact
```

大样本可显式启用实验性快速候选模式。它先用本项目自己的 canonical-k-mer
固定比例 sketch 建立内存倒排表，再对候选对调用与 exact 模式**完全相同**的
`estimate_with_index`；它不调用 skani，也不以 sketch 分数冒充 ANI/AF：

```bash
gtdb-ani-af triangle --list assemblies.list --out triangle.fast.tsv \
  --stats-out triangle.fast.stats.json \
  --threads 220 --profile medium --edge-mode sparse --min-af 15 \
  --triangle-mode sketch \
  --sketch-scale 64 \
  --triangle-min-shared 1 \
  --triangle-max-candidates 128 \
  --triangle-sketch-size 0 \
  --triangle-max-posting 1000 \
  --triangle-candidates-out triangle.candidates.tsv
```

- `sketch-scale` 决定固定比例采样密度，默认 64（约保留 1/64 canonical k-mer）；
  它同样必须进入 Lake 召回审计，较小值提高短重叠召回但增加索引体积；
- `triangle-min-shared` 是进入每个节点候选排名所需的最少共享 sketch hash，默认
  1；
- `triangle-max-candidates` 是每个节点按“共享 hash 降序、路径字典序”保留的
  top-K，默认 128；0 表示不截断；最终无向候选集合取两端 top-K 的并集；
- `triangle-sketch-size` 是每个基因组的 sketch 上限，默认 0（不截断）；
- `triangle-max-posting` 丢弃出现于过多输入基因组的高频 hash，默认 1000；0
  表示不丢弃；
- ANI/AF 行仍按原 list 的 `i`、`j` 顺序输出，候选线程完成顺序不影响结果。

只生成候选、不执行昂贵 ANI/AF 估计时，可省略 `--out` 并增加
`--triangle-candidates-only`。这个入口用于在 Lake 冻结正边上快速扫描不同的
top-K/共享数/高频 posting 参数，不会生成伪 ANI/AF：

```bash
gtdb-ani-af triangle --list assemblies.list --threads 220 \
  --triangle-mode sketch --triangle-candidates-only \
  --sketch-scale 64 \
  --triangle-min-shared 1 --triangle-max-candidates 128 \
  --triangle-sketch-size 0 --triangle-max-posting 1000 \
  --triangle-candidates-out triangle.candidates.tsv
```

快速候选只是高召回工程预筛，**不是零假阴性的数学保证**。冻结为生产参数前，
必须把 `triangle.candidates.tsv` 对照历史 exact/权威正边做召回审计。例如：

```bash
python3 tools/audit_triangle_candidate_recall.py \
  --candidates triangle.candidates.tsv \
  --truth skani_triangle_sparse.tsv \
  --manifest quality_pass_dna2bit_negative_sags.tsv \
  --skip-truth-outside-manifest \
  --min-ani 95 --min-max-af 15 --max-gc-diff 2 \
  --summary candidate_recall.json --missing candidate_misses.tsv \
  --require-recall 1.0
```

如果冻结的历史 triangle 覆盖的节点多于当前 manifest，必须显式给出
`--skip-truth-outside-manifest`；审计 JSON 会单独报告跳过行数，默认仍严格报错。
审计器只核对候选集合，不重算或修改 ANI/AF。对候选对产生的 ANI/AF 数值与
exact 路径逐字节相同；不在候选集合中的 pair 不会被计算，因此必须单独报告真实
Lake 召回率，不能把合成测试当成生产证明。

当前 Stage3B 固定 `--triangle-mode exact`，不会调用上述候选预筛接口。

search 使用索引冻结的 `k`、sketch scale、query sketch size 和 max posting；
显式传入不同值会报错。`--threads` 在 index 时并行读取/解析 reference，在 search
时并行处理 query，输出仍按 query-list 顺序确定。pair 和 search 都输出固定、
供现有解析器读取的五列表头：

```text
Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query
```

## 透明算法

### 1. FASTA 与方向

- 原生支持明文 FASTA 和 `.gz` FASTA；同一解压内容走完全相同的解析和计算路径。
- 支持 query/reference 的任意多个 contig；header 必须非空且唯一。
- 字母转大写，`U` 视作 `T`；非 A/C/G/T 的 IUPAC/N/其他字符统一视作不可
  比较的 `N`。
- N 会打断 seed，且不进入 ANI 或 AF 分母。AF 分母是该侧所有 contig 中可调用
  A/C/G/T 碱基数。
- mapper 明确分别扫描 query 正链与整条 contig 的反向互补链。

### 2. 数据库候选召回

每个完整 FASTA 枚举 canonical k-mer：取正向 2-bit k-mer 与其反向互补编码的
较小者，再经过固定公开的 SplitMix64 混合。`--sketch-scale S` 使用固定阈值
`floor(UINT64_MAX/S)`：reference 保存所有不超过阈值的去重 hash，query 使用
同一阈值、再把严格有序结果截到最多 `--query-sketch-size`（默认 1,000）。这是一
种自行实现的 FracMinHash 式固定比例抽样，不来自 skani。

同一阈值给出关键的子串性质：未突变 query k-mer 若被抽中，相同 hash 也必然在
完整 reference 的抽样中，不再受两侧独立 bottom-k 截断影响。候选分数叫
**query sketch hit fraction**（共享 hash 数 / 实际 query sketch 数）；它只是
排序证据，不被包装成 ANI、AF 或严格的统计 containment 估计。

建库把 `(hash, ref_id)` 批量排序后外部归并成 `POSTINGS.bin` 倒排表。search 对
每个 query hash 二分定位 postings，只给命中的 reference 累加计数；不会对每个
query 线性扫描全库的全部 reference sketch。只有排序后的前
`--max-candidates` 个 reference 才进入映射内核。出现于超过 `--max-posting`
个 reference（默认 1,000）的高频 hash 整组不写入倒排表，避免低复杂度/共有
序列扩大 postings 并污染排名。

默认 scale=64 时，均匀 hash 假设下每个 reference 约保存其不同 canonical
k-mer 的 1/64；5 Mb reference 约 78,125 条，199,923 个 5 Mb reference 在高频
过滤前的粗略上界约 156 亿条 posting，仅 32-bit ref id 就约 62.5 GB，尚未计
hash、offset、文件和内存开销。可见当前原型即使消除了逐 query 全库线性扫描，
仍未证明适合完整 GTDB 的建库资源。增大 scale 会减小索引但直接降低短 SAG
证据；不得未经真实召回对照自行调大。1 kb、高重复或更低 identity 的 SAG 仍
可能只有少量命中，必须测真实库上的 top-k recall。

索引还保留每个 reference 的严格有序 sketch、绝对路径、压缩文件大小/mtime、
解压后原始 FASTA 字节的 SHA-256、可调用碱基数和 sketch 计数。读取索引时会
校验 `COMPLETE.json` schema/status、三项 sketch 参数、reference 数、三个索引
文件 SHA-256、重复 reference、ref id 顺序、每条 sketch 严格递增、postings
hash/ref id 严格递增及逐 reference 计数闭合。

search 默认先用压缩文件的 size+mtime 做快速 metadata 校验；`--deep-verify`
会在候选 FASTA 本来就要读取和解析的同一次读入中，额外核验解压内容 SHA-256，
不会为哈希再预读一遍候选。需要审计运行时应启用 `--deep-verify`。

### 3. pair 映射、ANI 与双向 AF

1. reference 建立精确正向 k-mer occurrence 表；出现次数超过 `--max-occ`
   的重复 seed 丢弃。
2. query 每条 contig 的正链、反向互补链分别查表。命中按
   `(query_contig, reference_contig, strand, reference_pos-query_pos)` 投票。
3. 只保留至少 `--min-seed-hits` 票的 diagonal；在其重叠范围内按
   `--window` 切窗。
4. 窗内仅比较双方均为 A/C/G/T 的位置。至少 `--min-window` 个可比较碱基且
   identity 不低于 `--min-identity` 的窗成为候选映射。
5. 候选按 identity、可比较长度、seed 票数降序，坐标升序确定性排序。每个窗
   只收集 query 与 reference **双方都尚未占用**的 fresh 位置对；通过接受门后，
   计数和标记的是完全同一组 fresh 位置，保证一对一，不会一侧去重而另一侧重复。
6. `ANI = fresh matches / fresh aligned positions × 100`；由于同一 fresh pair
   同时覆盖两侧，`AF_query = fresh aligned / query callable bases`，
   `AF_ref = fresh aligned / reference callable bases`。

默认值：`k=15`、window 500、最小窗 100、最少 seed 票 3、最大 seed occurrence
64、最低窗 identity 70%。所有排名和 tie-break 都显式固定，因此同一输入和参数
输出确定。

## 重要限制

这是 ungapped diagonal/window 估计器，不是全局或局部动态规划 aligner：

- 密集 indel、结构重排、极短 contig、高重复区会造成覆盖低估或候选漏召回；
- 固定比例候选 seed 有限，远缘或覆盖极低 SAG 可能无法召回真实最近邻；
- ANI 是已接受、去重映射位置的 identity，不等于 skani 的 ANI 统计模型；
- AF 的 callable-base 分母及贪心去重口径也与 skani 不保证一致；
- `--min-identity` 是 mapper 接受窗的工程阈值，不应解释成物种边界。

五列接口兼容不代表数值兼容。当前 Stage 3B 不再做 reference 阳性/阴性判定；
所有通过质量门的 Dna2bit-negative SAG 都进入 exact SAG–SAG triangle。triangle
结果中的吸引边仍以 `ANI>=95` 为硬门并另加 GC 差护栏，AF 只进入边权。该接线
不能反向证明本引擎的 GTDB `index/search` 与 skani 数值等价。

## 索引提交与来源收据

index 只接受不存在的 `--out-dir`。它先在同一父目录创建全新的临时目录，生成
`REFS.tsv`、`SKETCHES.bin`、倒排 `POSTINGS.bin`，最后写
`COMPLETE.json`，再通过目录 rename 原子发布；失败会清理该次精确临时目录，
已有目标绝不覆盖。

`COMPLETE.json`（schema v3）绑定：schema/status、`k`、sketch scale、query sketch
size、max posting、reference/保留 posting 数，以及 reference manifest、taxonomy、species radii、release metadata、
`REFS.tsv`、`SKETCHES.bin`、`POSTINGS.bin` 的 SHA-256。FASTA 内容 SHA-256 是建库
时对**解压后的原始 FASTA 字节**冻结，所以同内容的 `.fna`/`.fna.gz` 计算一致，
同时仍能审计具体压缩源文件的 metadata。

## 编译与小型自测

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/gtdb-ani-af --self-test
```

内建自测使用确定性合成多 contig reference、正向 query 片段、反向互补 query
片段、约 1% substitutions 与 N gap，检查 ANI、双向一对一 AF、canonical
sketch、针对 3.1 Mb reference 的 1 kb/10 kb、95% identity 子集 top-1 候选召回、
gzip/明文等价、SHA-256 和重复运行确定性。集成小测也使用 3.1 Mb references，
实际建立含高频过滤的 postings 索引并执行多线程 search：

```bash
tests/gtdb_ani_af_synthetic.sh ./build/gtdb-ani-af
```

这些测试不会访问 GTDB 或启动大型任务，也不构成真实生产对照。
