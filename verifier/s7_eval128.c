/* Independent exact evaluator, __int128 Gram variant (copy of s7_eval.c with G stored as three int64 limb
 * matrices P, Q, R in the blob, G = 2^32 P + 2^16 Q + R recombined exactly in __int128 at load; magic 'SIG8').
 * Overflow: s7_prep128.py asserts 210*5040*M^2 + den7*5040*sum_b maxdiag(G_b) < 2^126, which bounds every partial
 * sum below (each per-block accumulator is a sum of at most #configurations <= 5040 Gram entries, |G| <= max diag).
 * Original header:  For every admissible 7-vertex r-graph H (r = 3 or 4; 35-bit colex mask):
 *   F(H)  = integer numerator of the objective (cosig3: 210 - sum_Q C(k_Q,2); sigma4: 24 e - sum_T c_T (c_T - 1)),
 *   Z(H)  = sum_blocks wb * sum over EVERY ordered root tuple theta with H[theta] == sigma_t (labelled, exact) and
 *           every ordered pair (U1, U2) of disjoint (m-s)-sets outside theta of G_t[flag(theta,U1), flag(theta,U2)],
 *           G_t = A_t^T A_t (integer), wb = 5040 / #configurations of the block,
 *   S(H)  = F * 5040 * M^2 - den7 * Z  =  den7 * 5040 * M^2 * (f(H) - sum_t <Q_t, M_t(H)>),   Q_t = G_t / M^2.
 * The ordered tuples theta are enumerated through a table built by s7_prep.py: for each sorted root set and its
 * labelled code c, the list of ALL orderings pi with code(pi) == canonical type mask (so every ordered theta with
 * H[theta] = sigma appears exactly once).  Unordered disjoint pairs {U1, U2} are counted as G[a][b] + G[b][a].
 * Arithmetic: __int128 accumulators, G entries int64 (|G| < 2^62 checked by the prep); no floating point.
 *
 * usage: s7_eval BLOB raw REPS.bin NTHREADS OUT.json            (all one-vertex extensions R | L << low)
 *        s7_eval BLOB list MASKS.bin NTHREADS OUT.json [EACH.txt] (every mask of the list; EACH.txt: per-mask values)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

typedef __int128 i128;
typedef int64_t i64;

static i64 *B; static size_t BN, BP;
static i64 rd(void) { if (BP >= BN) { fprintf(stderr, "blob underflow\n"); exit(2); } return B[BP++]; }
static i64 *rdarr(i64 *n) { *n = rd(); i64 *p = B + BP; BP += (size_t)*n; if (BP > BN) { fprintf(stderr, "blob arr\n"); exit(2); } return p; }

typedef struct { i64 nflags, *flagidx, goff; } Type;
typedef struct {
  i64 s, m, wb, nsig, nfree, nU, npairs, nrs, nperm, ntypes;
  i64 *pairs, *sigpos, *type_of, *moff, *mlist, *freepos;
  Type *T;
} Block;

static i64 R_, MODE, OBJ, DEN7, M2, NCOP, *COP, NOBJ, *OBJT, LOW, NLINK, NB;
static Block *BL; static i128 *G; static i64 GN;

static void load(const char *path) {
  FILE *f = fopen(path, "rb"); if (!f) { perror(path); exit(2); }
  fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  BN = (size_t)sz / 8; B = malloc((size_t)sz); if (fread(B, 8, BN, f) != BN) { fprintf(stderr, "read\n"); exit(2); } fclose(f);
  BP = 0;
  if (rd() != 0x53494738) { fprintf(stderr, "magic\n"); exit(2); }
  R_ = rd(); MODE = rd(); OBJ = rd(); DEN7 = rd(); M2 = rd();
  COP = rdarr(&NCOP); OBJT = rdarr(&NOBJ);
  LOW = rd(); NLINK = rd(); NB = rd();
  BL = calloc((size_t)NB, sizeof(Block));
  for (i64 b = 0; b < NB; b++) {
    Block *k = &BL[b]; i64 n;
    k->s = rd(); k->m = rd(); k->wb = rd(); k->nsig = rd(); k->nfree = rd(); k->nU = rd(); k->npairs = rd();
    k->nrs = rd(); k->nperm = rd(); k->ntypes = rd();
    k->pairs = rdarr(&n); if (n != 2 * k->npairs) exit(3);
    k->sigpos = rdarr(&n); if (n != k->nrs * k->nsig) exit(3);
    k->type_of = rdarr(&n); if (n != (1LL << k->nsig)) exit(3);
    k->moff = rdarr(&n); if (n != (1LL << k->nsig) + 1) exit(3);
    k->mlist = rdarr(&n); if (n != k->moff[1LL << k->nsig]) exit(3);
    k->freepos = rdarr(&n); if (n != k->nrs * k->nperm * k->nU * k->nfree) exit(3);
    if (k->nU > 64) exit(3);
    k->T = calloc((size_t)k->ntypes, sizeof(Type));
    for (i64 t = 0; t < k->ntypes; t++) {
      k->T[t].nflags = rd(); k->T[t].flagidx = rdarr(&n); if (n != (1LL << k->nfree)) exit(3); k->T[t].goff = rd();
    }
  }
  GN = rd();
  const i64 *Pl = B + BP, *Ql = B + BP + GN, *Rl = B + BP + 2 * GN; BP += 3 * (size_t)GN;
  if (BP != BN) { fprintf(stderr, "blob size mismatch %zu %zu\n", BP, BN); exit(2); }
  G = malloc((size_t)GN * sizeof(i128)); if (!G) { fprintf(stderr, "alloc G\n"); exit(2); }
  for (i64 i = 0; i < GN; i++) G[i] = (((i128)Pl[i]) << 32) + (((i128)Ql[i]) << 16) + (i128)Rl[i];
}

static int admissible(i64 h) {
  if (MODE == 0) { for (i64 i = 0; i < NCOP; i++) if ((h & COP[i]) == COP[i]) return 0; }
  else { for (i64 i = 0; i < NCOP; i++) if ((h & COP[i]) == 0) return 0; }
  return 1;
}

static i64 objF(i64 h) {
  if (OBJ == 0) { i64 s = 0; for (int q = 0; q < 35; q++) { i64 c = __builtin_popcountll(h & OBJT[q]); s += c * (c - 1) / 2; } return 210 - s; }
  i64 P = 0; for (int t = 0; t < 35; t++) { i64 c = __builtin_popcountll(h & OBJT[t]); P += c * (c - 1); }
  return 24 * __builtin_popcountll(h) - P;
}

/* returns 0 on success, 1 if a configuration has no type / flag (inadmissible sub-configuration) */
static int evalZ(i64 h, i128 *Zout) {
  i128 Z = 0;
  for (i64 b = 0; b < NB; b++) {
    const Block *k = &BL[b];
    i128 acc = 0;
    int a[64];
    for (i64 i = 0; i < k->nrs; i++) {
      i64 c = 0;
      const i64 *sp = k->sigpos + i * k->nsig;
      for (i64 j = 0; j < k->nsig; j++) c |= ((h >> sp[j]) & 1) << j;
      i64 t = k->type_of[c];
      if (t < 0) return 1;
      const Type *T = &k->T[t];
      const i128 *Gt = G + T->goff; i64 nf = T->nflags;
      for (i64 q = k->moff[c]; q < k->moff[c + 1]; q++) {
        i64 pi = k->mlist[q];
        const i64 *fp = k->freepos + ((i * k->nperm + pi) * k->nU) * k->nfree;
        for (i64 u = 0; u < k->nU; u++) {
          i64 fc = 0;
          for (i64 j = 0; j < k->nfree; j++) fc |= ((h >> fp[j]) & 1) << j;
          fp += k->nfree;
          i64 x = T->flagidx[fc];
          if (x < 0) return 1;
          a[u] = (int)x;
        }
        for (i64 p = 0; p < k->npairs; p++) {
          int x = a[k->pairs[2 * p]], y = a[k->pairs[2 * p + 1]];
          acc += Gt[(i64)x * nf + y] + Gt[(i64)y * nf + x];
        }
      }
    }
    Z += acc * k->wb;
  }
  *Zout = Z;
  return 0;
}

static void p128(FILE *f, i128 v) {
  char buf[64]; int n = 0, neg = v < 0; unsigned __int128 u = neg ? (unsigned __int128)(-v) : (unsigned __int128)v;
  if (u == 0) buf[n++] = '0';
  while (u) { buf[n++] = '0' + (int)(u % 10); u /= 10; }
  if (neg) fputc('-', f);
  while (n) fputc(buf[--n], f);
}

#define KLOW 256
#define NFMAX 1300
typedef struct {
  i64 nadm, nrej, nbad;
  i128 maxZ[NFMAX]; i64 argZ[NFMAX], cnt[NFMAX];
  i128 low[KLOW]; i64 lowm[KLOW], lowc[KLOW]; int nlow, wmax;
} Acc;

static i64 FMIN, NF;

static void acc_init(Acc *A) {
  memset(A, 0, sizeof(*A));
  for (int g = 0; g < NFMAX; g++) A->maxZ[g] = -((i128)1 << 126);
}

static void low_insert(Acc *A, i128 S, i64 h, i64 c) {
  if (A->nlow == KLOW && S > A->low[A->wmax]) return;
  for (int i = 0; i < A->nlow; i++) if (A->low[i] == S) { A->lowc[i] += c; if (h < A->lowm[i]) A->lowm[i] = h; return; }
  int w;
  if (A->nlow < KLOW) w = A->nlow++;
  else w = A->wmax;
  A->low[w] = S; A->lowm[w] = h; A->lowc[w] = c;
  A->wmax = 0; for (int i = 1; i < A->nlow; i++) if (A->low[i] > A->low[A->wmax]) A->wmax = i;
}

static void acc_add(Acc *A, i64 h, i64 F, i128 Z) {
  i64 g = F - FMIN; if (g < 0 || g >= NF) { fprintf(stderr, "F range\n"); exit(4); }
  A->cnt[g]++;
  if (Z > A->maxZ[g] || (Z == A->maxZ[g] && h < A->argZ[g])) { A->maxZ[g] = Z; A->argZ[g] = h; }
  i128 S = (i128)F * 5040 * M2 - (i128)DEN7 * Z;
  low_insert(A, S, h, 1);
}

static void acc_merge(Acc *D, const Acc *A) {
  D->nadm += A->nadm; D->nrej += A->nrej; D->nbad += A->nbad;
  for (int g = 0; g < NFMAX; g++) {
    D->cnt[g] += A->cnt[g];
    if (A->cnt[g] && (A->maxZ[g] > D->maxZ[g] || (A->maxZ[g] == D->maxZ[g] && A->argZ[g] < D->argZ[g]))) {
      D->maxZ[g] = A->maxZ[g]; D->argZ[g] = A->argZ[g];
    }
  }
  for (int i = 0; i < A->nlow; i++) low_insert(D, A->low[i], A->lowm[i], A->lowc[i]);
}

int main(int argc, char **argv) {
  if (argc < 6) { fprintf(stderr, "usage\n"); return 1; }
  load(argv[1]);
  const char *mode = argv[2];
  int nth = atoi(argv[4]);
  omp_set_num_threads(nth);
  if (OBJ == 0) { FMIN = 0; NF = 211; } else { FMIN = -420; NF = 1261; }
  FILE *f = fopen(argv[3], "rb"); if (!f) { perror(argv[3]); return 2; }
  fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
  i64 nm = sz / 8; i64 *ms = malloc((size_t)sz); if (fread(ms, 8, (size_t)nm, f) != (size_t)nm) return 2; fclose(f);
  Acc *D = malloc(sizeof(Acc)); acc_init(D);
  i64 *ext = calloc((size_t)nm, sizeof(i64));
  FILE *each = (argc > 6) ? fopen(argv[6], "w") : NULL;
  double t0 = omp_get_wtime();
  if (!strcmp(mode, "raw")) {
    i64 nl = 1LL << NLINK, chunk = 1LL << 12, nper = nl / chunk, ntask = nm * nper;
    #pragma omp parallel
    {
      Acc *A = malloc(sizeof(Acc)); acc_init(A);
      #pragma omp for schedule(dynamic, 1)
      for (i64 task = 0; task < ntask; task++) {
        i64 ri = task / nper, l0 = (task % nper) * chunk, e = 0;
        for (i64 L = l0; L < l0 + chunk; L++) {
          i64 h = ms[ri] | (L << LOW);
          if (!admissible(h)) { A->nrej++; continue; }
          i128 Z; if (evalZ(h, &Z)) { A->nbad++; continue; }
          A->nadm++; e++;
          acc_add(A, h, objF(h), Z);
        }
        #pragma omp atomic
        ext[ri] += e;
      }
      #pragma omp critical
      acc_merge(D, A);
      free(A);
    }
  } else {
    i128 *Zs = calloc((size_t)nm, sizeof(i128)); int *st = calloc((size_t)nm, sizeof(int));
    #pragma omp parallel
    {
      Acc *A = malloc(sizeof(Acc)); acc_init(A);
      #pragma omp for schedule(dynamic, 1024)
      for (i64 i = 0; i < nm; i++) {
        i64 h = ms[i];
        if (!admissible(h)) { A->nrej++; st[i] = 2; continue; }
        i128 Z; if (evalZ(h, &Z)) { A->nbad++; st[i] = 1; continue; }
        A->nadm++; Zs[i] = Z;
        acc_add(A, h, objF(h), Z);
      }
      #pragma omp critical
      acc_merge(D, A);
      free(A);
    }
    if (each) {
      for (i64 i = 0; i < nm; i++) {
        fprintf(each, "%lld %d %lld ", (long long)ms[i], st[i], (long long)(st[i] ? 0 : objF(ms[i])));
        p128(each, st[i] ? 0 : Zs[i]); fputc(' ', each);
        p128(each, st[i] ? 0 : (i128)objF(ms[i]) * 5040 * M2 - (i128)DEN7 * Zs[i]); fputc('\n', each);
      }
      fclose(each);
    }
  }
  double t1 = omp_get_wtime();
  /* sort lowest */
  for (int i = 0; i < D->nlow; i++) for (int j = i + 1; j < D->nlow; j++)
    if (D->low[j] < D->low[i] || (D->low[j] == D->low[i] && D->lowm[j] < D->lowm[i])) {
      i128 t = D->low[i]; D->low[i] = D->low[j]; D->low[j] = t; i64 u = D->lowm[i]; D->lowm[i] = D->lowm[j]; D->lowm[j] = u;
      u = D->lowc[i]; D->lowc[i] = D->lowc[j]; D->lowc[j] = u; }
  FILE *o = fopen(argv[5], "w");
  fprintf(o, "{\"mode\": \"%s\", \"inputs\": %lld, \"admissible\": %lld, \"rejected\": %lld, \"bad\": %lld, \"seconds\": %.1f, \"threads\": %d,\n",
          mode, (long long)nm, (long long)D->nadm, (long long)D->nrej, (long long)D->nbad, t1 - t0, nth);
  fprintf(o, " \"S_den\": "); p128(o, (i128)DEN7 * 5040 * M2); fprintf(o, ",\n \"lowest\": [");
  for (int i = 0; i < D->nlow; i++) { fprintf(o, "%s[\"", i ? ", " : ""); p128(o, D->low[i]); fprintf(o, "\", %lld, %lld]", (long long)D->lowm[i], (long long)D->lowc[i]); }
  fprintf(o, "],\n \"groups\": [");
  int first = 1;
  for (int g = 0; g < NF; g++) if (D->cnt[g]) {
    fprintf(o, "%s[%lld, %lld, \"", first ? "" : ", ", (long long)(g + FMIN), (long long)D->cnt[g]); p128(o, D->maxZ[g]);
    fprintf(o, "\", %lld]", (long long)D->argZ[g]); first = 0;
  }
  fprintf(o, "]");
  if (!strcmp(mode, "raw")) { fprintf(o, ",\n \"ext\": ["); for (i64 i = 0; i < nm; i++) fprintf(o, "%s%lld", i ? ", " : "", (long long)ext[i]); fprintf(o, "]"); }
  fprintf(o, "}\n");
  fclose(o);
  fprintf(stderr, "done %s: admissible %lld rejected %lld bad %lld in %.1fs\n", mode, (long long)D->nadm, (long long)D->nrej, (long long)D->nbad, t1 - t0);
  return 0;
}
