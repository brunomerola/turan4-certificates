/* Independent exact evaluator (own code, adapted from the design of rv_eval.c; new objective,
 * multi-limb exact arithmetic, unordered-pair symmetry, per-F groups).  Reads the blob written by rs_prep.py.
 *
 * For a 4-graph H on {0..6} (own encoding: bit i = i-th 4-subset in LEXICOGRAPHIC order):
 *   F(H) = 24 e(H) - P(H),  P(H) = sum over 3-sets T of c_T (c_T - 1)   (c_T = #{x not in T : T+x in H}),
 *   Zb(H) = sum over ORDERED injective s-tuples theta with H[theta] == sigma_k (labelled, exact test)
 *           sum over ORDERED pairs (U1, U2) of disjoint (m-s)-subsets of V - theta   G_k[flag(theta,U1), flag(theta,U2)]
 *         = 2 * sum over UNORDERED pairs {U1, U2}  (G_k symmetric, checked exactly in rs_prep.py),
 *   Zsum(H) = sum_b (5040 / conf_b) Zb(H),  conf_b = 7!/(7-s)! C(7-s, m-s) C(7-m, m-s),
 * and with G_k = D Q_k:  slack(H) = f(H) - sum_b <Q_b, M_b(H)> = F/420 - Zsum/(5040 D).  Reported exactly:
 *   S(H) = 5040 D (slack(H) - b) = 12 D F - 5040 D b - Zsum      (b = target of the blob).
 * Arithmetic: every G entry is stored as NL balanced limbs of 45 bits (|limb| <= 2^44); per graph each limb is
 * accumulated in int64 (bound sum_b 5040 * 2^44 < 2^60, checked in rs_prep.py), then S is normalised with exact
 * carries to digits d_0..d_{NL-2} in [0, 2^45) and a signed top digit; comparisons are lexicographic (exact).
 *
 * Modes (single thread; split work over processes with LO HI):
 *   rs_eval BLOB list  MASKS.bin LO HI [TIGHT.bin]   MASKS.bin = uint64 COLEX masks (producer class list payload)
 *   rs_eval BLOB raw   REPS.bin  LO HI [TIGHT.bin]   REPS.bin = uint64 own masks of 6-vertex reps (vertex 6 unused);
 *                                                    tasks LO..HI-1 of nreps*1024 (rep, high 10 link bits); all 2^20
 *                                                    links of vertex 6 per rep (completeness by heredity)
 *   rs_eval BLOB each  MASKS.bin                     prints S for every graph (cross-checks)
 * TIGHT.bin receives the own masks of all graphs with S == 0.
 *
 * -DSYM variant (binary rs_eval4s): blob from rs_prep.py --sym-out, whose matrices are Gsym_k = sum_{h in Aut(sigma_k)}
 * G_k[pi_h ., pi_h .].  For every UNORDERED root set S it takes the FIRST ordering theta0 of S with H[theta0] equal to
 * a key's sigma (the types of one block are pairwise non-isomorphic, asserted in rs_prep.py), and adds
 * Gsym[flag(theta0,U1), flag(theta0,U2)] once: the valid orderings of S are exactly theta0 o h, h in Aut(sigma), and
 * flag(theta0 o h, U) = pi_{h^-1} flag(theta0, U).  Same S(H) as the brute-force binary; checked graph by graph.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <sys/mman.h>

#ifndef NLC
#define NLC 4
#endif
#define LB 45
typedef long long ll;
typedef __int128 i128;

#define NV 7
#define NE 35
static int E7[NE][4];
static int colex_of_own[NE], own_of_colex[NE];
static int TRI[35][4];   /* for each 3-set T (lex), own indices of T+x, x not in T */

typedef struct {
    ll s, m, nsig, nE, conf, nkeys;
    ll *V, *T;
    int nth, nU, npairs, w;
    int *rootpos, *pos;
    int pairs[64][2];
    int active;
    int nS, sfact, *Sth;       /* SYM: nS root sets, each with sfact orderings (theta indices) */
} Block;

static int nblocks, nkeys;
static Block B[8];
static ll *Qk[64];
static ll Qn[64];
static int Qnz[64];
static ll c12L[NLC], TtL[NLC];
static uint64_t pmask[21];

static int comb(int n, int k) { if (k < 0 || k > n) return 0; int r = 1; for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1); return r; }

static int own_index(int a, int b, int c, int d) {
    int v[4] = {a, b, c, d};
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (v[j] < v[i]) { int t = v[i]; v[i] = v[j]; v[j] = t; }
    for (int i = 0; i < NE; i++) if (E7[i][0] == v[0] && E7[i][1] == v[1] && E7[i][2] == v[2] && E7[i][3] == v[3]) return i;
    fprintf(stderr, "bad edge\n"); exit(1);
}

static void init_edges(void) {
    int n = 0;
    for (int a = 0; a < 7; a++) for (int b = a + 1; b < 7; b++) for (int c = b + 1; c < 7; c++) for (int d = c + 1; d < 7; d++) {
        E7[n][0] = a; E7[n][1] = b; E7[n][2] = c; E7[n][3] = d;
        colex_of_own[n] = comb(a, 1) + comb(b, 2) + comb(c, 3) + comb(d, 4);
        n++;
    }
    for (int i = 0; i < NE; i++) own_of_colex[colex_of_own[i]] = i;
    int t = 0;
    for (int a = 0; a < 7; a++) for (int b = a + 1; b < 7; b++) for (int c = b + 1; c < 7; c++) {
        int q = 0;
        for (int x = 0; x < 7; x++) if (x != a && x != b && x != c) TRI[t][q++] = own_index(a, b, c, x);
        t++;
    }
    int p = 0;
    for (int S = 0; S < 128; S++) if (__builtin_popcount(S) == 5) {
        uint64_t mk = 0;
        for (int i = 0; i < NE; i++) { int in = 1; for (int j = 0; j < 4; j++) if (!(S >> E7[i][j] & 1)) in = 0; if (in) mk |= 1ULL << i; }
        pmask[p++] = mk;
    }
    if (p != 21) { fprintf(stderr, "psets\n"); exit(1); }
}

static uint64_t colex_to_own(uint64_t h) { uint64_t o = 0; for (int i = 0; i < NE; i++) if ((h >> i) & 1) o |= 1ULL << own_of_colex[i]; return o; }
static uint64_t own_to_colex(uint64_t h) { uint64_t o = 0; for (int i = 0; i < NE; i++) if ((h >> i) & 1) o |= 1ULL << colex_of_own[i]; return o; }
static void rd(void *p, size_t sz, size_t n, FILE *f) { if (fread(p, sz, n, f) != n) { fprintf(stderr, "short read\n"); exit(1); } }

static void load_blob(const char *fn) {
    FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
    ll hdr[5]; rd(hdr, 8, 5, f);
#ifdef SYM
    const ll MAGIC = 0x53594D31;
#else
    const ll MAGIC = 0x52375347;
#endif
    if (hdr[0] != MAGIC || hdr[3] != NLC || hdr[4] != LB) { fprintf(stderr, "blob header / NLC mismatch (NL=%lld)\n", hdr[3]); exit(1); }
    nblocks = (int)hdr[1]; nkeys = (int)hdr[2];
    rd(c12L, 8, NLC, f); rd(TtL, 8, NLC, f);
    for (int b = 0; b < nblocks; b++) {
        Block *X = &B[b];
        ll h[6]; rd(h, 8, 6, f);
        X->s = h[0]; X->m = h[1]; X->nsig = h[2]; X->nE = h[3]; X->conf = h[4]; X->nkeys = h[5];
        ll *order = malloc(sizeof(ll) * X->nE * 4); rd(order, 8, X->nE * 4, f);
        X->V = malloc(sizeof(ll) << X->nsig); rd(X->V, 8, 1ULL << X->nsig, f);
        X->T = malloc(sizeof(ll) << X->nE); rd(X->T, 8, 1ULL << X->nE, f);
        int s = (int)X->s, m = (int)X->m, k = m - s;
        if (2 * m - s > NV || 5040 % X->conf) { fprintf(stderr, "bad block\n"); exit(1); }
        int nth = 1; for (int i = 0; i < s; i++) nth *= (NV - i);
        X->nth = nth; X->nU = comb(NV - s, k);
        X->rootpos = malloc(sizeof(int) * nth * (X->nsig ? X->nsig : 1));
        X->pos = malloc(sizeof(int) * nth * X->nU * X->nE);
        int th[8], t = 0, tot = 1;
        for (int i = 0; i < s; i++) tot *= NV;
        int *usedm = malloc(sizeof(int) * nth);
        for (int code = 0; code < tot; code++) {
            int c = code, used = 0, ok = 1;
            for (int i = 0; i < s; i++) { th[i] = c % NV; c /= NV; if (used >> th[i] & 1) ok = 0; used |= 1 << th[i]; }
            if (!ok) continue;
            int rest[7], nr = 0; for (int v = 0; v < NV; v++) if (!(used >> v & 1)) rest[nr++] = v;
            for (int j = 0; j < X->nsig; j++)
                X->rootpos[t * X->nsig + j] = own_index(th[order[4*j]], th[order[4*j+1]], th[order[4*j+2]], th[order[4*j+3]]);
            int u = 0;
            for (int um = 0; um < (1 << nr); um++) {
                if (__builtin_popcount(um) != k) continue;
                int lab[8]; for (int i = 0; i < s; i++) lab[i] = th[i];
                int q = s; for (int i = 0; i < nr; i++) if (um >> i & 1) lab[q++] = rest[i];
                for (int j = 0; j < X->nE; j++)
                    X->pos[(t * X->nU + u) * X->nE + j] = own_index(lab[order[4*j]], lab[order[4*j+1]], lab[order[4*j+2]], lab[order[4*j+3]]);
                u++;
            }
            if (u != X->nU) { fprintf(stderr, "nU\n"); exit(1); }
            usedm[t] = used;
            t++;
        }
        if (t != nth) { fprintf(stderr, "nth\n"); exit(1); }
        /* group the ordered root tuples by their unordered root set */
        X->nS = comb(NV, s); X->sfact = 1; for (int i = 2; i <= s; i++) X->sfact *= i;
        X->Sth = malloc(sizeof(int) * X->nS * X->sfact);
        { int si = 0;
          for (int Sm = 0; Sm < (1 << NV); Sm++) if (__builtin_popcount(Sm) == s) {
              int c2 = 0;
              for (int q = 0; q < nth; q++) if (usedm[q] == Sm) X->Sth[si * X->sfact + c2++] = q;
              if (c2 != X->sfact) { fprintf(stderr, "Sth\n"); exit(1); }
              si++;
          }
          if (si != X->nS) { fprintf(stderr, "nS\n"); exit(1); } }
        free(usedm);
        int nr = NV - s, ul[64], nu = 0;
        for (int um = 0; um < (1 << nr); um++) if (__builtin_popcount(um) == k) ul[nu++] = um;
        X->npairs = 0;
        int nordered = 0;
        for (int a = 0; a < nu; a++) for (int c = 0; c < nu; c++) if (a != c && !(ul[a] & ul[c])) {
            nordered++;
            if (a < c) { X->pairs[X->npairs][0] = a; X->pairs[X->npairs][1] = c; X->npairs++; }
        }
        if ((ll)nth * nordered != X->conf || nordered != 2 * X->npairs) { fprintf(stderr, "conf mismatch\n"); exit(1); }
        X->w = (int)(2 * 5040 / X->conf);     /* ordered pairs = 2 x unordered; weight 5040/conf */
        if ((ll)nth * X->npairs > 8192) { fprintf(stderr, "entry buffer too small\n"); exit(1); }
        free(order);
    }
    for (int i = 0; i < nkeys; i++) {
        ll h[3]; rd(h, 8, 3, f);
        int k = (int)h[0]; Qn[k] = h[1]; Qnz[k] = (int)h[2];
        { size_t bytes = sizeof(ll) * h[1] * h[1] * NLC, al = 1 << 21;
          Qk[k] = aligned_alloc(al, ((bytes + al - 1) / al) * al);
          madvise(Qk[k], ((bytes + al - 1) / al) * al, MADV_HUGEPAGE); }   /* fewer TLB misses (performance only) */
        rd(Qk[k], 8, h[1] * h[1] * NLC, f);
    }
    char extra; if (fread(&extra, 1, 1, f) == 1) { fprintf(stderr, "trailing bytes in blob\n"); exit(1); }
    fclose(f);
    for (int b = 0; b < nblocks; b++) {
        B[b].active = 0;
        for (ll r = 0; r < (1LL << B[b].nsig); r++) if (B[b].V[r] >= 0 && Qnz[B[b].V[r]]) B[b].active = 1;
    }
}

static inline int admissible(uint64_t h) { for (int i = 0; i < 21; i++) if (!(h & pmask[i])) return 0; return 1; }

static inline int objF(uint64_t h) {
    int P = 0;
    for (int t = 0; t < 35; t++) {
        int c = (int)((h >> TRI[t][0]) & 1) + (int)((h >> TRI[t][1]) & 1) + (int)((h >> TRI[t][2]) & 1) + (int)((h >> TRI[t][3]) & 1);
        P += c * (c - 1);
    }
    return 24 * __builtin_popcountll(h) - P;
}

/* S normalised: S[0..NLC-2] in [0, 2^45), S[NLC-1] signed */
#ifdef TWOPHASE
static const ll *ebuf[8192];    /* >= max over blocks of (#ordered root tuples) x (#unordered pairs) = 2520 */
#endif
static inline void eval_S(uint64_t h, int F, ll *S) {
    ll Zs[NLC];
    for (int j = 0; j < NLC; j++) Zs[j] = 0;
    int fl[64];
    for (int b = 0; b < nblocks; b++) {
        Block *X = &B[b];
        if (!X->active) continue;
        ll Zb[NLC];
        for (int j = 0; j < NLC; j++) Zb[j] = 0;
#ifdef TWOPHASE
        int ne = 0;     /* phase 1 collects the entry addresses (prefetched), phase 2 sums them */
#endif
        int nsig = (int)X->nsig, nE = (int)X->nE, nU = X->nU;
#ifdef SYM
        for (int si = 0; si < X->nS; si++) {
          for (int ti = 0; ti < X->sfact; ti++) {
            int t = X->Sth[si * X->sfact + ti];
            ll r = 0;
            const int *rp = X->rootpos + t * nsig;
            for (int j = 0; j < nsig; j++) r |= (ll)((h >> rp[j]) & 1) << j;
            ll k = X->V[r];
            if (k < 0) continue;
            if (!Qnz[k]) break;
#else
        for (int t = 0; t < X->nth; t++) {
            ll r = 0;
            const int *rp = X->rootpos + t * nsig;
            for (int j = 0; j < nsig; j++) r |= (ll)((h >> rp[j]) & 1) << j;
            ll k = X->V[r];
            if (k < 0 || !Qnz[k]) continue;
#endif
            const int *pp = X->pos + (ll)t * nU * nE;
            for (int u = 0; u < nU; u++) {
                ll x = r;
                const int *q = pp + u * nE;
                for (int j = nsig; j < nE; j++) x |= (ll)((h >> q[j]) & 1) << j;
                fl[u] = (int)X->T[x];
            }
            const ll *Q = Qk[k]; ll n = Qn[k];
            for (int q = 0; q < X->npairs; q++) {
                int a = fl[X->pairs[q][0]], c = fl[X->pairs[q][1]];
                if (a >= 0 && c >= 0) {
                    const ll *e = Q + ((ll)a * n + c) * NLC;
#ifdef TWOPHASE
                    __builtin_prefetch(e);
                    ebuf[ne++] = e;
#else
                    for (int j = 0; j < NLC; j++) Zb[j] += e[j];
#endif
                }
            }
#ifdef SYM
            break;
          }
#endif
        }
#ifdef TWOPHASE
        for (int i = 0; i < ne; i++) { const ll *e = ebuf[i]; for (int j = 0; j < NLC; j++) Zb[j] += e[j]; }
#endif
        for (int j = 0; j < NLC; j++) Zs[j] += (ll)X->w * Zb[j];
    }
    ll carry = 0;
    for (int j = 0; j < NLC; j++) {
        ll v = c12L[j] * F - TtL[j] - Zs[j] + carry;
        if (j < NLC - 1) {
            ll d = v & ((1LL << LB) - 1);
            carry = (v - d) >> LB;      /* exact: v - d is a multiple of 2^45 */
            S[j] = d;
        } else S[j] = v;
    }
}

static inline int cmpS(const ll *A, const ll *Bv) {
    for (int j = NLC - 1; j >= 0; j--) if (A[j] != Bv[j]) return A[j] < Bv[j] ? -1 : 1;
    return 0;
}
static inline int isZero(const ll *A) { for (int j = 0; j < NLC; j++) if (A[j]) return 0; return 1; }
static inline int isNeg(const ll *A) { return A[NLC - 1] < 0; }

#define KLOW 256
typedef struct {
    ll cnt[841], zc[841]; ll mn[841][NLC]; uint64_t am[841];
    ll nadm, nrej, nzero, nneg;
    ll gmin[NLC]; uint64_t hmin;
    ll pmin[NLC]; uint64_t hp; int havep;      /* smallest POSITIVE S */
    ll low[KLOW][NLC]; uint64_t lowh[KLOW]; int nlow;
    FILE *tight;
} Acc;

static void acc_init(Acc *A) { memset(A, 0, sizeof(Acc)); }

static void lowins(Acc *A, const ll *S, uint64_t h) {
    if (A->nlow < KLOW) { memcpy(A->low[A->nlow], S, sizeof(ll) * NLC); A->lowh[A->nlow] = h; A->nlow++; return; }
    int w = 0; for (int i = 1; i < KLOW; i++) if (cmpS(A->low[i], A->low[w]) > 0) w = i;
    if (cmpS(S, A->low[w]) < 0) { memcpy(A->low[w], S, sizeof(ll) * NLC); A->lowh[w] = h; }
}

static int keep_low = 1;
static inline void process(Acc *A, uint64_t h) {
    if (!admissible(h)) { A->nrej++; return; }
    int F = objF(h);
    if (F < 0 || F > 840) { fprintf(stderr, "F range\n"); exit(1); }
    ll S[NLC];
    eval_S(h, F, S);
    A->nadm++;
    if (A->cnt[F] == 0 || cmpS(S, A->mn[F]) < 0) { memcpy(A->mn[F], S, sizeof(S)); A->am[F] = h; }
    A->cnt[F]++;
    if (A->nadm == 1 || cmpS(S, A->gmin) < 0) { memcpy(A->gmin, S, sizeof(S)); A->hmin = h; }
    if (isZero(S)) { A->nzero++; A->zc[F]++; if (A->tight) fwrite(&h, 8, 1, A->tight); }
    else if (isNeg(S)) { A->nneg++; }
    else if (!A->havep || cmpS(S, A->pmin) < 0) { memcpy(A->pmin, S, sizeof(S)); A->hp = h; A->havep = 1; }
    if (keep_low) lowins(A, S, h);
}

static void prS(const char *tag, const ll *S, uint64_t h) {
    printf("%s", tag);
    for (int j = 0; j < NLC; j++) printf(" %lld", S[j]);
    printf(" | own %llu colex %llu\n", (unsigned long long)h, (unsigned long long)own_to_colex(h));
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
    init_edges();
    load_blob(argv[1]);
    printf("# NLC=%d nblocks=%d nkeys=%d\n", NLC, nblocks, nkeys);
    for (int b = 0; b < nblocks; b++) printf("# block s=%lld m=%lld conf=%lld nth=%d nU=%d unordered_pairs=%d w=%d active=%d\n", B[b].s, B[b].m, B[b].conf, B[b].nth, B[b].nU, B[b].npairs, B[b].w, B[b].active);
    fflush(stdout);
    Acc *A = malloc(sizeof(Acc)); acc_init(A);
    clock_t t0 = clock();
    FILE *f = fopen(argv[3], "rb"); if (!f) { perror(argv[3]); return 1; }
    fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
    ll n = sz / 8; uint64_t *H = malloc(sz); rd(H, 8, n, f); fclose(f);
    if (!strcmp(argv[2], "each")) {
        for (ll i = 0; i < n; i++) {
            uint64_t h = colex_to_own(H[i]);
            if (!admissible(h)) { printf("EACH %llu inadmissible\n", (unsigned long long)H[i]); continue; }
            ll S[NLC]; int F = objF(h); eval_S(h, F, S);
            printf("EACH colex %llu F %d S", (unsigned long long)H[i], F);
            for (int j = 0; j < NLC; j++) printf(" %lld", S[j]);
            printf("\n");
        }
        return 0;
    }
    ll lo = atoll(argv[4]), hi = atoll(argv[5]);
    if (argc > 6) A->tight = fopen(argv[6], "wb");
    if (!strcmp(argv[2], "list")) {
        if (hi > n) hi = n;
        for (ll i = lo; i < hi; i++) process(A, colex_to_own(H[i]));
    } else if (!strcmp(argv[2], "raw")) {
        keep_low = 0;
        int lpos[20], j = 0;
        for (int a = 0; a < 6; a++) for (int b = a + 1; b < 6; b++) for (int c = b + 1; c < 6; c++) lpos[j++] = own_index(a, b, c, 6);
        uint64_t lw[1024], hg[1024];
        for (int x = 0; x < 1024; x++) { lw[x] = hg[x] = 0; for (int i = 0; i < 10; i++) if (x >> i & 1) { lw[x] |= 1ULL << lpos[i]; hg[x] |= 1ULL << lpos[10 + i]; } }
        for (ll r = 0; r < n; r++) for (int i = 0; i < NE; i++) if ((H[r] >> i & 1) && E7[i][3] == 6) { fprintf(stderr, "rep uses vertex 6\n"); return 1; }
        ll ntask = n * 1024;
        if (hi > ntask) hi = ntask;
        for (ll task = lo; task < hi; task++) {
            uint64_t base = H[task / 1024] | hg[task % 1024];
            for (int x = 0; x < 1024; x++) process(A, base | lw[x]);
            if ((task - lo) % 4096 == 4095) { fprintf(stderr, "# progress %lld/%lld %.0fs\n", task - lo + 1, hi - lo, (double)(clock() - t0) / CLOCKS_PER_SEC); }
        }
    } else { fprintf(stderr, "mode\n"); return 1; }
    if (A->tight) fclose(A->tight);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("RANGE %lld %lld\n", lo, hi);
    printf("TOTAL admissible=%lld rejected=%lld zero=%lld negative=%lld cpu=%.1fs\n", A->nadm, A->nrej, A->nzero, A->nneg, secs);
    prS("MIN", A->gmin, A->hmin);
    if (A->havep) prS("NEXT", A->pmin, A->hp);
    for (int F = 0; F <= 840; F++) if (A->cnt[F]) {
        printf("G %d %lld %lld", F, A->cnt[F], A->zc[F]);
        prS("", A->mn[F], A->am[F]);
    }
    for (int i = 0; i < A->nlow; i++) prS("LOW", A->low[i], A->lowh[i]);
    return 0;
}
