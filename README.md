# Microsags v0.1（DNA2bit-SAG Tractor original）

把以下流程统一封装成一个 C++17 命令行程序：

1. 主入口读取解压后的文件内容自动识别输入，不依赖扩展名：
   - 成对 FASTQ：从 `fastp` 质控开始；
   - FASTA contigs：跳过 `fastp` 和 SPAdes，直接从长度门开始；
   - 单端 FASTQ、R1/R2 混合格式、损坏或截断文件均 fail-closed；
2. reads 路径执行 `SPAdes --sc --careful`；contigs 路径保持原序列不变（未压缩文件建立链接，gzip 输入仅解压为稳定 FASTA 视图）；
3. 两条路径统一统计每个SAG组装序列总长度；`<1000 bp` 硬排除，`>=1000 bp` 才继续；
4. 使用老师原版 `/home/data/shared/software/Dna2bit/dna2bit`，读取该SAG质控后的配对FASTQ给每个合格SAG分类；
5. 有标签的 SAG 按 species 标签聚合，进入已有 C++ `cpp-subass-full`；
6. 无标签 SAG 写入 `03B_unclassified_pending.tsv`；独立二进制
   `sag-stage3b-tractor` 可继续执行严格的无标签 SAG 聚合流程（见
   [README_STAGE3B.md](README_STAGE3B.md)）。主程序与 3B 可分别部署，避免改变
   已冻结的 1–3A 行为。

这里“C++封装”是指主调度、清单检查、分组、FASTA合并、收据和3B接口均为C++；fastp、SPAdes、老师的dna2bit及Flye科学内核仍作为固定依赖调用，不冒充重新实现这些成熟算法。

## 冻结的 dna2bit 规则

- 原程序：`/home/data/shared/software/Dna2bit/dna2bit`
- GTDB R232 bit库：`/home/data/shared/software/Dna2bit/GTDB232/GTDB`
- 与该bit库一一对应的taxonomy：`/home/data/shared/software/Dna2bit/GTDB232/genome_taxonomy_1.csv`
  （199,923行；服务器现存文件SHA-256：`b10a9a80d42455f7fdc419afcbf231fb06119ee5de4346d5cb5bf8ba9749b78b`）
- `k=17`
- `bit_len=55296`（必须与现有参考bit一致）
- `hash_type=0`（wyhash）
- reads 路径：严格沿用老师历史流程，对质控后 R1/R2 执行 `sketch -p`
- contigs 路径：老师原版 `sketch` 直接读取 FASTA，不加 `-p`；其余参数完全相同
- `min_ratio=0.01`（老师原版search接受门）
- 不使用Cellbit57、ALC、Top16、GC硬门或训练模型。
- 对老师二进制的目录参数强制保留末尾 `/`；原程序内部直接拼接文件名，缺少它会段错误或找错路径。
- 老师二进制还存在长路径段错误：封装器在每个工作目录建立不改数据的短名符号链接，并从该目录用相对短路径调用；输出再由 manifest/收据绑定回原 SAG。该兼容层不改变任何序列或 sketch 参数。
- SAG长度硬门：全部contig碱基总和 `<1000 bp` 排除，恰好1000 bp保留。

参考库小文件加载的严格兼容优化候选见
[docs/DNA2BIT_PACKED_SEARCH.md](docs/DNA2BIT_PACKED_SEARCH.md)。主流程默认仍为
`--dna-search-engine teacher`，不会改变旧行为；只有完成全量逐字节验收后，才可显式使用
`--dna-search-engine packed --dna-packed-db INDEX_DIR`。生产 packed index 在构建时把正确
taxonomy的完整 accession 集合与SHA绑定进去，搜索时错配旧taxonomy会立即失败，而不会
静默产生大量空标签。

## 输入

TSV 可以混合两种 SAG 输入；程序逐行读取内容自动分流。表头可有可无。

成对 reads（三列）：

```text
SAG_ID<TAB>/absolute/R1.fastq.gz<TAB>/absolute/R2.fastq.gz
```

已有 contigs（两列）：

```text
SAG_ID<TAB>/absolute/assembly.fasta[.gz]
```

检测依据是解压后的 FASTA/FASTQ 结构，不是 `.fa`、`.fq`、`.gz` 等文件名。相对路径按 manifest 所在目录解析。

## 编译和运行

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/dna2bit-sag-pipeline --manifest SAGs.tsv --out output --threads 48 --memory-gb 256
```

可先加 `--dry-run` 查看命令；意外中断后加 `--resume`，只复用带PASS收据的阶段。

## 输出

- `01_assembly/<SAG>/`：fastp报告、SPAdes结果和阶段收据
- `INPUT_AUDIT.tsv`：逐 SAG 记录检测类型、原始输入、自动起点、跳过步骤、统一 assembly、总 bp 与长度门结果
- `01_assembly/excluded_lt1000bp.tsv`：未达到1000 bp、不会进入后续步骤的SAG
- `STAGE3B_FASTA_STATS.tsv`：在长度门同一次FASTA扫描中取得的 total/max-contig/GC，
  供新跑的CheckM2结果按 `sag_id` 合并；不伪造 contamination，也不复用历史值
- `02_dna2bit/bits/`：每个SAG的老师原版bit
- `02_dna2bit/labels.tsv`：被原版门接受的标签
- `03A_subassemble/<species>/run/assembly.fasta`：有标签路径最终bin
- 3A 显式使用 `cpp-subass --no-overlap-policy passthrough`。它只在
  `flye-modules assemble` 返回 0、draft 为真实 0 字节文件且本次专用日志明确报告
  `No overlaps found`/`Assembled 0 disjointigs` 时，把已验证的聚合输入逐字节保留为
  `assembly.fasta`；单 SAG 与多 SAG 组因此都不会只因“彼此无 overlap”丢失 bin。
  每个这类 bin 都有 `NO_OVERLAP_PASSTHROUGH.PASS.json` 记录命令、触发证据、
  输入/输出 SHA-256、记录/碱基/GC/长度统计以及 ID/记录集合闭包。其余错误仍直接失败；
  非空 draft 始终进入原 Flye repeat/contigger/polishing 路径。
- `03B_unclassified_pending.tsv`：留给未来3B插件
- `TIMING.tsv` / `TIMING.json`：steady clock记录的 preflight+length、sketch、search、
  3A subassemble 和 total 墙钟时间
- `COMPLETE.json`：全流程闭包

## 当前限制

主程序本身仍止于 3B pending 接口；3B 由独立可执行文件实现。当前不做任务并发
调度、资源门或容器化。species 字符串直接来自老师原版 dna2bit 结果，不加入
后续魔改规则。

仓库另含实验性的独立 C++17 ANI/AF 估计器 `gtdb-ani-af`，详见
[README_GTD_ANI_AF.md](README_GTD_ANI_AF.md)。它不调用 skani，也未经生产对照，
Stage 3B 仅在显式 `--allow-experimental-ani-engine` 时接入。当前 3B 不搜索 GTDB
reference、不接收 `--gtdb-taxonomy`/`--ani-db`，而是让所有通过
`max_contig>=1000`、CheckM2 `contamination<5` 的 Dna2bit-negative SAG 进入
exact SAG–SAG triangle；pairwise 吸引边仍要求 `ANI>=95` 并通过 GC 差护栏。
