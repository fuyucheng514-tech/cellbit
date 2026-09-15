# Dna2bit packed search（老师版严格兼容候选）

## 目的与当前边界

老师部署版 `dna2bit search` 每次先用 `opendir/readdir` 枚举 GTDB bit
目录，再逐个 `fopen("rb")` 读取 199,923 个 6,912-byte 小文件。这个候选只把
参考库的物理布局改成一个连续文件；距离、接受门和输出语义均按部署二进制逆向
结果保留。它目前是独立可执行程序，尚未替换生产主流程。

## 已冻结的可观察语义

- 参考和 query 都使用原始 `opendir/readdir` 顺序，不排序。
- 每个 bit 由 864 个本机 `uint64_t` 组成；距离为逐词
  `popcount(query XOR reference)` 之和。
- 只有严格更小的距离才更新 winner，所以距离相等时目录顺序更早者胜出。
- 老师版的 `second` 并不是真正的全局第二小距离：每次出现严格的新纪录最小时，
  仅把旧 `best` 赋给 `second`。最终 `second` 是最终 winner 出现前的上一条纪录。
- margin 用单精度 `float` 计算 `(second-best)/best`；只有严格
  `margin > 0.01` 才接受，恰好 `0.01` 拒绝。
- 文件名用老师二进制内置正则 `.*(GC[FA]_\d+\.\d+)_` 提取 accession。
- 拒绝项不写 CSV；接受项为 `query完整路径,accession,taxonomy\n`。老师的 CRLF taxonomy
  读取方式会保留 `\r`，本实现原样保留，因此输出可逐字节一致。

## 编译与小型对照测试

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j --target dna2bit-pack-builder dna2bit-packed-search
python3 tests/test_dna2bit_packed_teacher_equivalence.py \
  --teacher /home/data/shared/software/Dna2bit/dna2bit \
  --packed "$PWD/build/dna2bit-packed"
```

测试只创建 2–3 个参考的临时目录，不访问正式 GTDB 库。它覆盖错误但必须兼容的
previous-record `second`、严格阈值、严格 tie、zero distance 和单参考边界，并要求
老师二进制与 packed search 的 CSV 逐字节相等。

对照命令刻意使用老师程序的短参数 `-d/-s/-t/-m/-n/-o`。部署二进制虽然在 help
中列出了长参数，但实测全长参数组合会段错误；这也是现有主包装器使用短参数的原因。

## 一次性构建及搜索（全库暂未执行）

```bash
dna2bit-pack-builder \
  --database /path/to/GTDB/ \
  --index /path/to/GTDB232.k17.packed \
  --taxonomy /path/to/GTDB232/genome_taxonomy_1.csv

dna2bit-packed-search \
  --index /path/to/GTDB232.k17.packed \
  --search /path/to/query_bits/ \
  --min-ratio 0.01 --nthreads 48 \
  --result-file search_result.csv
```

目录参数末尾 `/` 是兼容契约的一部分。builder 按当次 raw `readdir` 顺序冻结
`references.pack`、`references.tsv`、两者各自的 SHA-256 与 PASS receipt，然后原子发布目录。
search 默认在 mmap/计算前同时核验 manifest 和 pack 的 SHA-256；任何清单或数据位翻转
都会 fail-closed，除非仅为受控诊断显式传入 `--no-verify-pack-sha`。

生产构建必须绑定 taxonomy。显式传 `--taxonomy` 最清楚；若省略，builder 只会自动寻找
bit目录同级的 `genome_taxonomy_1.csv`。builder 会要求 taxonomy 与参考 manifest 的 accession
集合严格一一对应，再把原始字节流及其 SHA-256 绑定进索引；search 未传 `--tax` 时默认
使用这份内置表，若显式传入另一份表则必须与绑定 SHA 完全相同。这样 199,923 条
GTDB232 reference 不可能再误配 113,104 行的旧 taxonomy；检测到恰好199,923条参考却
仍找不到taxonomy时会直接失败。无同级表且小于全库规模的未绑定 developer index 仅保留
给微型兼容性测试，搜索时必须显式传 `--tax`。
正式 199,923 行 pack 的精确大小为 1,381,867,776 bytes（约 1.287 GiB），应作为
外置数据库，不计入小于 1 GiB 的程序包。

## 正式替换前仍须完成

1. 在服务器空闲时构建一次全库 packed index，冻结源目录 inventory 与顺序。
2. 先抽 64–256 条 Lake query，与老师版 `-m 0.01` 逐字节比对。
3. 用 1、8、220 线程验证 packed 输出完全一致。
4. 最后对全部 8,785 条 Lake query 做逐字节验收；通过前不得切换生产入口。
