#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <zlib.h>

// Isolated pair-only research candidate. It is deliberately not wired into Stage 3B.
// No learned correction and no skani code/data are used here.
namespace fs = std::filesystem;

struct Config {
  int k = 15;
  int seed_scale = 4;
  int max_occ = 32;
  int min_anchors = 3;
  int max_predecessors = 256;
  int max_chain_gap = 12000;
  int max_diag_drift = 160;
  int align_band = 96;
  int min_aligned = 100;
  double min_identity = 70.0;
  fs::path query, reference, out;
};

struct Contig { std::string name, seq; };
struct Genome {
  fs::path path;
  std::vector<Contig> contigs;
  std::uint64_t callable_bases = 0;
};
struct Occ { std::uint32_t contig = 0, pos = 0; };
struct Anchor { std::uint32_t q = 0, r = 0; };
struct GroupKey {
  std::uint32_t qc = 0, rc = 0;
  bool reverse = false;
  bool operator<(const GroupKey &o) const {
    return std::tie(qc, rc, reverse) < std::tie(o.qc, o.rc, o.reverse);
  }
};
struct Operation {
  std::int64_t q = -1, r = -1; // Oriented query and forward reference coordinates.
  char code = 0;               // M, X, I(query base), D(reference base).
};
struct ChainAlignment {
  GroupKey key;
  int anchor_count = 0, chain_score = 0;
  std::uint32_t q0 = 0, q1 = 0, r0 = 0, r1 = 0;
  std::vector<Operation> ops;
  std::uint64_t matches = 0, denominator = 0, q_covered = 0, r_covered = 0;
  double identity = 0.0;
};
struct Estimate {
  std::uint64_t matches = 0, mismatches = 0, insertions = 0, deletions = 0;
  std::uint64_t q_covered = 0, r_covered = 0;
  double ani = 0.0, af_query = 0.0, af_ref = 0.0;
  int accepted_chains = 0;
};

static int base_code(char c) {
  switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return -1; }
}
static bool callable(char c) { return base_code(c) >= 0; }
static char complement(char c) {
  switch (c) { case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; }
}
static std::string reverse_complement(const std::string &s) {
  std::string r(s.size(), 'N');
  for (std::size_t i = 0; i < s.size(); ++i) r[s.size() - 1 - i] = complement(s[i]);
  return r;
}
static std::uint64_t mix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}
static std::string trim(std::string s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
  return s;
}

static std::string read_text(const fs::path &p) {
  if (p.extension() == ".gz") {
    gzFile z = gzopen(p.string().c_str(), "rb");
    if (!z) throw std::runtime_error("cannot open gzip FASTA: " + p.string());
    std::string out;
    char buf[1 << 16];
    int n = 0;
    while ((n = gzread(z, buf, sizeof(buf))) > 0) out.append(buf, static_cast<std::size_t>(n));
    int err = Z_OK;
    const char *msg = gzerror(z, &err);
    if (err != Z_OK && err != Z_STREAM_END) {
      std::string why = msg ? msg : "gzip read error";
      gzclose(z);
      throw std::runtime_error(why);
    }
    if (gzclose(z) != Z_OK) throw std::runtime_error("gzip checksum/close failure: " + p.string());
    return out;
  }
  std::ifstream in(p, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open FASTA: " + p.string());
  std::ostringstream out;
  out << in.rdbuf();
  return out.str();
}

static Genome read_fasta(const fs::path &p) {
  if (!fs::is_regular_file(p)) throw std::runtime_error("FASTA missing: " + p.string());
  Genome g;
  g.path = fs::absolute(p).lexically_normal();
  std::istringstream in(read_text(p));
  std::string line;
  Contig *cur = nullptr;
  std::set<std::string> names;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line.front() == '>') {
      std::string name = trim(line.substr(1));
      if (name.empty() || !names.insert(name).second) throw std::runtime_error("empty/duplicate FASTA header");
      g.contigs.push_back({name, {}});
      cur = &g.contigs.back();
      continue;
    }
    if (!cur && !trim(line).empty()) throw std::runtime_error("sequence before FASTA header");
    if (!cur) continue;
    for (unsigned char x : line) if (!std::isspace(x)) {
      char c = static_cast<char>(std::toupper(x));
      if (c == 'U') c = 'T';
      if (!callable(c)) c = 'N'; else ++g.callable_bases;
      cur->seq.push_back(c);
    }
  }
  g.contigs.erase(std::remove_if(g.contigs.begin(), g.contigs.end(), [](const Contig &c) { return c.seq.empty(); }), g.contigs.end());
  if (g.contigs.empty() || g.callable_bases == 0) throw std::runtime_error("FASTA has no callable bases");
  return g;
}

template <class Fn>
static void each_kmer(const std::string &s, int k, Fn fn) {
  const std::uint64_t mask = (1ULL << (2 * k)) - 1ULL;
  std::uint64_t code = 0;
  int valid = 0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    int b = base_code(s[i]);
    if (b < 0) { code = 0; valid = 0; continue; }
    code = ((code << 2) | static_cast<std::uint64_t>(b)) & mask;
    if (++valid >= k) fn(static_cast<std::uint32_t>(i - static_cast<std::size_t>(k) + 1), code);
  }
}

using SeedIndex = std::unordered_map<std::uint64_t, std::vector<Occ>>;
static SeedIndex build_reference_seeds(const Genome &reference, const Config &cfg) {
  SeedIndex index;
  for (std::size_t ci = 0; ci < reference.contigs.size(); ++ci) {
    each_kmer(reference.contigs[ci].seq, cfg.k, [&](std::uint32_t pos, std::uint64_t code) {
      if (mix64(code) % static_cast<std::uint64_t>(cfg.seed_scale) != 0) return;
      auto &v = index[code];
      if (static_cast<int>(v.size()) <= cfg.max_occ) v.push_back({static_cast<std::uint32_t>(ci), pos});
    });
  }
  for (auto it = index.begin(); it != index.end();) {
    if (static_cast<int>(it->second.size()) > cfg.max_occ) it = index.erase(it); else ++it;
  }
  return index;
}

static std::map<GroupKey, std::vector<Anchor>> collect_anchors(const Genome &query, const SeedIndex &index, const Config &cfg) {
  std::map<GroupKey, std::vector<Anchor>> groups;
  for (std::size_t qi = 0; qi < query.contigs.size(); ++qi) {
    for (int strand = 0; strand < 2; ++strand) {
      const std::string oriented = strand ? reverse_complement(query.contigs[qi].seq) : query.contigs[qi].seq;
      each_kmer(oriented, cfg.k, [&](std::uint32_t qpos, std::uint64_t code) {
        if (mix64(code) % static_cast<std::uint64_t>(cfg.seed_scale) != 0) return;
        auto found = index.find(code);
        if (found == index.end()) return;
        for (const Occ &o : found->second) groups[{static_cast<std::uint32_t>(qi), o.contig, strand != 0}].push_back({qpos, o.pos});
      });
    }
  }
  return groups;
}

// Sparse deterministic local collinear chaining. This explicitly allows the
// diagonal to drift, so substitutions and short indels do not fragment a hit.
static std::vector<std::vector<Anchor>> chain_group(std::vector<Anchor> a, const Config &cfg) {
  std::sort(a.begin(), a.end(), [](const Anchor &x, const Anchor &y) { return std::tie(x.q, x.r) < std::tie(y.q, y.r); });
  a.erase(std::unique(a.begin(), a.end(), [](const Anchor &x, const Anchor &y) { return x.q == y.q && x.r == y.r; }), a.end());
  const std::size_t n = a.size();
  std::vector<int> score(n, cfg.k * 10), pred(n, -1);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t begin = i > static_cast<std::size_t>(cfg.max_predecessors) ? i - static_cast<std::size_t>(cfg.max_predecessors) : 0;
    for (std::size_t j = i; j-- > begin;) {
      if (a[i].q <= a[j].q || a[i].r <= a[j].r) continue;
      const int dq = static_cast<int>(a[i].q - a[j].q), dr = static_cast<int>(a[i].r - a[j].r);
      if (dq > cfg.max_chain_gap || dr > cfg.max_chain_gap) continue;
      const int drift = std::abs(dq - dr);
      const int allowed = cfg.max_diag_drift + std::max(dq, dr) / 20;
      if (drift > allowed) continue;
      const int nonoverlap = std::min({cfg.k, dq, dr});
      const int gap_penalty = drift * 2 + static_cast<int>(std::log2(static_cast<double>(std::max(dq, dr)) + 1.0));
      const int candidate = score[j] + nonoverlap * 10 - gap_penalty;
      if (candidate > score[i] || (candidate == score[i] && static_cast<int>(j) < pred[i])) {
        score[i] = candidate;
        pred[i] = static_cast<int>(j);
      }
    }
  }
  std::vector<std::size_t> endpoints(n);
  for (std::size_t i = 0; i < n; ++i) endpoints[i] = i;
  std::sort(endpoints.begin(), endpoints.end(), [&](std::size_t x, std::size_t y) {
    if (score[x] != score[y]) return score[x] > score[y];
    return std::tie(a[x].q, a[x].r) < std::tie(a[y].q, a[y].r);
  });
  std::vector<unsigned char> used(n, 0);
  std::vector<std::vector<Anchor>> chains;
  for (std::size_t end : endpoints) {
    std::vector<std::size_t> path;
    for (int at = static_cast<int>(end); at >= 0; at = pred[static_cast<std::size_t>(at)]) path.push_back(static_cast<std::size_t>(at));
    std::reverse(path.begin(), path.end());
    int fresh = 0;
    for (std::size_t x : path) fresh += !used[x];
    if (fresh < cfg.min_anchors || fresh * 2 < static_cast<int>(path.size())) continue;
    std::vector<Anchor> chain;
    for (std::size_t x : path) {
      if (!used[x]) chain.push_back(a[x]);
    }
    // Forced exact anchors must not overlap one another.
    std::vector<Anchor> spaced;
    for (const Anchor &x : chain) {
      if (spaced.empty() || (x.q >= spaced.back().q + static_cast<std::uint32_t>(cfg.k) && x.r >= spaced.back().r + static_cast<std::uint32_t>(cfg.k))) spaced.push_back(x);
    }
    if (static_cast<int>(spaced.size()) < cfg.min_anchors) continue;
    if (spaced.back().q + cfg.k - spaced.front().q < static_cast<std::uint32_t>(cfg.min_aligned)) continue;
    for (std::size_t x : path) used[x] = 1;
    chains.push_back(std::move(spaced));
  }
  return chains;
}

struct TraceRow {
  int lo = 0, hi = -1;
  std::vector<int> score;
  std::vector<unsigned char> move; // 0 diag, 1 up(I), 2 left(D), 3 unreachable.
  int get(int j) const {
    if (j < lo || j > hi) return std::numeric_limits<int>::min() / 8;
    return score[static_cast<std::size_t>(j - lo)];
  }
};

// Global alignment of the interval between two exact anchors. The DP is
// restricted to a band around the straight line connecting those anchors.
static std::vector<Operation> align_gap(const std::string &q, std::uint32_t qb, std::uint32_t qe,
                                        const std::string &r, std::uint32_t rb, std::uint32_t re,
                                        int requested_band) {
  const int n = static_cast<int>(qe - qb), m = static_cast<int>(re - rb);
  std::vector<Operation> out;
  if (n == 0) { for (int j = 0; j < m; ++j) out.push_back({-1, static_cast<std::int64_t>(rb + j), 'D'}); return out; }
  if (m == 0) { for (int i = 0; i < n; ++i) out.push_back({static_cast<std::int64_t>(qb + i), -1, 'I'}); return out; }
  const int band = std::max(requested_band, std::abs(n - m) + 8);
  const int neg = std::numeric_limits<int>::min() / 8;
  std::vector<TraceRow> rows(static_cast<std::size_t>(n + 1));
  for (int i = 0; i <= n; ++i) {
    const int center = static_cast<int>(std::llround(static_cast<double>(i) * m / n));
    TraceRow &row = rows[static_cast<std::size_t>(i)];
    row.lo = std::max(0, center - band);
    row.hi = std::min(m, center + band);
    row.score.assign(static_cast<std::size_t>(row.hi - row.lo + 1), neg);
    row.move.assign(row.score.size(), 3);
    for (int j = row.lo; j <= row.hi; ++j) {
      const std::size_t x = static_cast<std::size_t>(j - row.lo);
      if (i == 0 && j == 0) { row.score[x] = 0; continue; }
      int best = neg;
      unsigned char move = 3;
      if (i > 0 && j > 0) {
        int prev = rows[static_cast<std::size_t>(i - 1)].get(j - 1);
        if (prev > neg / 2) {
          const char qc = q[static_cast<std::size_t>(qb) + i - 1], rc = r[static_cast<std::size_t>(rb) + j - 1];
          best = prev + ((callable(qc) && callable(rc) && qc == rc) ? 2 : -3);
          move = 0;
        }
      }
      if (i > 0) {
        int prev = rows[static_cast<std::size_t>(i - 1)].get(j);
        if (prev > neg / 2 && prev - 4 > best) { best = prev - 4; move = 1; }
      }
      if (j > row.lo && row.score[x - 1] > neg / 2 && row.score[x - 1] - 4 > best) {
        best = row.score[x - 1] - 4;
        move = 2;
      }
      row.score[x] = best;
      row.move[x] = move;
    }
  }
  if (rows.back().get(m) <= neg / 2) throw std::runtime_error("banded alignment endpoint unreachable");
  int i = n, j = m;
  while (i > 0 || j > 0) {
    const TraceRow &row = rows[static_cast<std::size_t>(i)];
    if (j < row.lo || j > row.hi) throw std::runtime_error("banded traceback escaped band");
    const unsigned char move = row.move[static_cast<std::size_t>(j - row.lo)];
    if (move == 0) {
      --i; --j;
      char code = q[static_cast<std::size_t>(qb) + i] == r[static_cast<std::size_t>(rb) + j] ? 'M' : 'X';
      out.push_back({static_cast<std::int64_t>(qb + i), static_cast<std::int64_t>(rb + j), code});
    } else if (move == 1) {
      --i;
      out.push_back({static_cast<std::int64_t>(qb + i), -1, 'I'});
    } else if (move == 2) {
      --j;
      out.push_back({-1, static_cast<std::int64_t>(rb + j), 'D'});
    } else throw std::runtime_error("invalid banded traceback");
  }
  std::reverse(out.begin(), out.end());
  return out;
}

static ChainAlignment align_chain(const GroupKey &key, const std::vector<Anchor> &chain,
                                  const Genome &query, const Genome &reference, const Config &cfg) {
  const std::string q = key.reverse ? reverse_complement(query.contigs[key.qc].seq) : query.contigs[key.qc].seq;
  const std::string &r = reference.contigs[key.rc].seq;
  ChainAlignment a;
  a.key = key;
  a.anchor_count = static_cast<int>(chain.size());
  a.q0 = chain.front().q; a.r0 = chain.front().r;
  a.q1 = chain.back().q + static_cast<std::uint32_t>(cfg.k);
  a.r1 = chain.back().r + static_cast<std::uint32_t>(cfg.k);
  auto add_seed = [&](const Anchor &x) {
    for (int z = 0; z < cfg.k; ++z) {
      const std::uint32_t qp = x.q + static_cast<std::uint32_t>(z), rp = x.r + static_cast<std::uint32_t>(z);
      a.ops.push_back({qp, rp, q[qp] == r[rp] ? 'M' : 'X'});
    }
  };
  add_seed(chain.front());
  for (std::size_t i = 1; i < chain.size(); ++i) {
    const Anchor &prev = chain[i - 1], &next = chain[i];
    const std::uint32_t qb = prev.q + static_cast<std::uint32_t>(cfg.k), rb = prev.r + static_cast<std::uint32_t>(cfg.k);
    std::vector<Operation> middle = align_gap(q, qb, next.q, r, rb, next.r, cfg.align_band);
    a.ops.insert(a.ops.end(), middle.begin(), middle.end());
    add_seed(next);
  }
  for (const Operation &op : a.ops) {
    if (op.q >= 0 && callable(q[static_cast<std::size_t>(op.q)])) ++a.q_covered;
    if (op.r >= 0 && callable(r[static_cast<std::size_t>(op.r)])) ++a.r_covered;
    if (op.q >= 0 && op.r >= 0) {
      const char qc = q[static_cast<std::size_t>(op.q)], rc = r[static_cast<std::size_t>(op.r)];
      if (!callable(qc) || !callable(rc)) continue;
      ++a.denominator;
      if (qc == rc) ++a.matches;
    } else if (op.q >= 0) {
      if (callable(q[static_cast<std::size_t>(op.q)])) ++a.denominator;
    } else if (op.r >= 0 && callable(r[static_cast<std::size_t>(op.r)])) ++a.denominator;
  }
  a.identity = a.denominator ? 100.0 * static_cast<double>(a.matches) / a.denominator : 0.0;
  a.chain_score = a.anchor_count * cfg.k * 10;
  return a;
}

static Estimate estimate(const Genome &query, const Genome &reference, const Config &cfg) {
  SeedIndex index = build_reference_seeds(reference, cfg);
  auto groups = collect_anchors(query, index, cfg);
  std::vector<ChainAlignment> alignments;
  for (auto &[key, anchors] : groups) {
    for (const auto &chain : chain_group(std::move(anchors), cfg)) {
      ChainAlignment a = align_chain(key, chain, query, reference, cfg);
      if (a.denominator >= static_cast<std::uint64_t>(cfg.min_aligned) && a.identity >= cfg.min_identity) alignments.push_back(std::move(a));
    }
  }
  std::sort(alignments.begin(), alignments.end(), [](const ChainAlignment &a, const ChainAlignment &b) {
    if (a.identity != b.identity) return a.identity > b.identity;
    if (a.denominator != b.denominator) return a.denominator > b.denominator;
    if (a.anchor_count != b.anchor_count) return a.anchor_count > b.anchor_count;
    return std::tie(a.key.qc, a.key.rc, a.key.reverse, a.q0, a.r0) < std::tie(b.key.qc, b.key.rc, b.key.reverse, b.q0, b.r0);
  });
  std::vector<std::vector<unsigned char>> qcov, rcov;
  for (const Contig &c : query.contigs) qcov.emplace_back(c.seq.size(), 0);
  for (const Contig &c : reference.contigs) rcov.emplace_back(c.seq.size(), 0);
  Estimate result;
  for (const ChainAlignment &a : alignments) {
    const std::size_t qlen = query.contigs[a.key.qc].seq.size();
    std::uint64_t possible = 0, fresh = 0;
    for (const Operation &op : a.ops) {
      std::int64_t oq = op.q;
      if (oq >= 0 && a.key.reverse) oq = static_cast<std::int64_t>(qlen - 1) - oq;
      bool qok = oq >= 0 && callable(query.contigs[a.key.qc].seq[static_cast<std::size_t>(oq)]);
      bool rok = op.r >= 0 && callable(reference.contigs[a.key.rc].seq[static_cast<std::size_t>(op.r)]);
      if (qok || rok) {
        ++possible;
        if ((!qok || !qcov[a.key.qc][static_cast<std::size_t>(oq)]) && (!rok || !rcov[a.key.rc][static_cast<std::size_t>(op.r)])) ++fresh;
      }
    }
    if (fresh < static_cast<std::uint64_t>(cfg.min_aligned) || fresh * 2 < possible) continue;
    ++result.accepted_chains;
    for (const Operation &op : a.ops) {
      std::int64_t oq = op.q;
      if (oq >= 0 && a.key.reverse) oq = static_cast<std::int64_t>(qlen - 1) - oq;
      bool qok = oq >= 0 && callable(query.contigs[a.key.qc].seq[static_cast<std::size_t>(oq)]);
      bool rok = op.r >= 0 && callable(reference.contigs[a.key.rc].seq[static_cast<std::size_t>(op.r)]);
      if ((qok && qcov[a.key.qc][static_cast<std::size_t>(oq)]) || (rok && rcov[a.key.rc][static_cast<std::size_t>(op.r)])) continue;
      if (qok) qcov[a.key.qc][static_cast<std::size_t>(oq)] = 1;
      if (rok) rcov[a.key.rc][static_cast<std::size_t>(op.r)] = 1;
      if (op.q >= 0 && op.r >= 0) {
        const char original_qc = query.contigs[a.key.qc].seq[static_cast<std::size_t>(oq)];
        const char qc = a.key.reverse ? complement(original_qc) : original_qc;
        const char rc = reference.contigs[a.key.rc].seq[static_cast<std::size_t>(op.r)];
        if (!qok || !rok) continue;
        if (qc == rc) ++result.matches; else ++result.mismatches;
      } else if (qok) ++result.insertions;
      else if (rok) ++result.deletions;
    }
  }
  for (std::size_t i = 0; i < qcov.size(); ++i) for (std::size_t j = 0; j < qcov[i].size(); ++j) result.q_covered += qcov[i][j] && callable(query.contigs[i].seq[j]);
  for (std::size_t i = 0; i < rcov.size(); ++i) for (std::size_t j = 0; j < rcov[i].size(); ++j) result.r_covered += rcov[i][j] && callable(reference.contigs[i].seq[j]);
  const std::uint64_t denominator = result.matches + result.mismatches + result.insertions + result.deletions;
  result.ani = denominator ? 100.0 * static_cast<double>(result.matches) / denominator : 0.0;
  result.af_query = 100.0 * static_cast<double>(result.q_covered) / query.callable_bases;
  result.af_ref = 100.0 * static_cast<double>(result.r_covered) / reference.callable_bases;
  return result;
}

static void write_result(std::ostream &out, const Genome &reference, const Genome &query, const Estimate &e) {
  out << "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tM\tX\tI\tD\tAccepted_chains\n";
  out << reference.path.string() << '\t' << query.path.string() << '\t'
      << std::fixed << std::setprecision(6) << e.ani << '\t' << e.af_ref << '\t' << e.af_query << '\t'
      << e.matches << '\t' << e.mismatches << '\t' << e.insertions << '\t' << e.deletions << '\t' << e.accepted_chains << '\n';
}

static std::string mutate_with_indels(std::string s) {
  for (std::size_t i = 41; i < s.size(); i += 97) s[i] = "ACGT"[(base_code(s[i]) + 1) % 4];
  s.insert(900, "ACGTTGCA");
  s.erase(1700, 11);
  return s;
}
static std::string deterministic_dna(std::size_t n) {
  std::string s(n, 'A');
  std::uint64_t x = 0x123456789abcdef0ULL;
  for (char &c : s) { x ^= x << 7; x ^= x >> 9; x ^= x << 8; c = "ACGT"[x & 3ULL]; }
  return s;
}
static void self_test() {
  Config cfg;
  cfg.k = 13;
  cfg.seed_scale = 2;
  cfg.min_aligned = 80;
  std::string reference_seq = deterministic_dna(7000);
  std::string query_seq = mutate_with_indels(reference_seq.substr(1200, 3200));
  Genome reference{{"synthetic-reference"}, {{"r", reference_seq}}, static_cast<std::uint64_t>(reference_seq.size())};
  Genome query{{"synthetic-query"}, {{"q", query_seq}}, static_cast<std::uint64_t>(query_seq.size())};
  Estimate a = estimate(query, reference, cfg), b = estimate(query, reference, cfg);
  if (std::tie(a.matches, a.mismatches, a.insertions, a.deletions, a.q_covered, a.r_covered) !=
      std::tie(b.matches, b.mismatches, b.insertions, b.deletions, b.q_covered, b.r_covered)) throw std::runtime_error("non-deterministic result");
  if (a.ani < 97.0 || a.af_query < 90.0 || a.af_ref < 40.0 || a.insertions + a.deletions == 0) throw std::runtime_error("gapped alignment self-test failed");
  std::cout << "PASS gapped pair self-test ANI=" << a.ani << " AFq=" << a.af_query << " AFr=" << a.af_ref
            << " I+D=" << (a.insertions + a.deletions) << '\n';
}

static void help() {
  std::cout << "gtdb-ani-af-gapped-candidate (isolated pair-only prototype)\n"
            << "  pair --query Q.fa[.gz] --reference R.fa[.gz] [--out result.tsv]\n"
            << "  --k 15 --seed-scale 4 --max-occ 32 --min-anchors 3\n"
            << "  --max-chain-gap 12000 --align-band 96 --min-aligned 100 --min-identity 70\n"
            << "  --self-test\n";
}

int main(int argc, char **argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--self-test") { self_test(); return 0; }
    if (argc < 2 || std::string(argv[1]) != "pair") { help(); return argc < 2 ? 0 : 2; }
    Config cfg;
    for (int i = 2; i < argc; ++i) {
      std::string arg = argv[i];
      auto value = [&]() -> std::string { if (++i >= argc) throw std::runtime_error("missing value for " + arg); return argv[i]; };
      if (arg == "--query") cfg.query = value();
      else if (arg == "--reference") cfg.reference = value();
      else if (arg == "--out") cfg.out = value();
      else if (arg == "--k") cfg.k = std::stoi(value());
      else if (arg == "--seed-scale") cfg.seed_scale = std::stoi(value());
      else if (arg == "--max-occ") cfg.max_occ = std::stoi(value());
      else if (arg == "--min-anchors") cfg.min_anchors = std::stoi(value());
      else if (arg == "--max-chain-gap") cfg.max_chain_gap = std::stoi(value());
      else if (arg == "--align-band") cfg.align_band = std::stoi(value());
      else if (arg == "--min-aligned") cfg.min_aligned = std::stoi(value());
      else if (arg == "--min-identity") cfg.min_identity = std::stod(value());
      else if (arg == "--help" || arg == "-h") { help(); return 0; }
      else throw std::runtime_error("unknown option: " + arg);
    }
    if (cfg.query.empty() || cfg.reference.empty()) throw std::runtime_error("pair requires --query and --reference");
    if (cfg.k < 5 || cfg.k > 31 || cfg.seed_scale < 1 || cfg.max_occ < 1 || cfg.min_anchors < 2 ||
        cfg.max_chain_gap < cfg.k || cfg.align_band < 8 || cfg.min_aligned < cfg.k || cfg.min_identity < 0 || cfg.min_identity > 100) throw std::runtime_error("invalid numeric option");
    Genome query = read_fasta(cfg.query), reference = read_fasta(cfg.reference);
    Estimate result = estimate(query, reference, cfg);
    if (cfg.out.empty()) write_result(std::cout, reference, query, result);
    else {
      if (fs::exists(cfg.out)) throw std::runtime_error("refusing to overwrite output: " + cfg.out.string());
      std::ofstream out(cfg.out);
      if (!out) throw std::runtime_error("cannot create output: " + cfg.out.string());
      write_result(out, reference, query, result);
      if (!out) throw std::runtime_error("failed writing output: " + cfg.out.string());
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "ERROR: " << e.what() << '\n';
    return 1;
  }
}
