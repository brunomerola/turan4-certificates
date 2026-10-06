/* Independent evaluator (own code; reads the blob written by rv_prep.py).
 *
 * For a 4-graph H on vertex set {0..6} (own encoding: bit i = i-th 4-subset in lexicographic order) computes, from
 * the definition of the flag-algebra averaging,
 *     Zb(H) = sum over ORDERED injective root tuples theta (s-tuples) with H[theta] == sigma_k (labelled, exactly)
 *             sum over ORDERED pairs (U1, U2) of disjoint (m-s)-subsets of V - theta
 *             Qint_k[ flag(theta, U1), flag(theta, U2) ]          (0 if a flag is not in key k's list)
 *     Zsum(H) = sum_b (5040 / conf_b) Zb(H),   conf_b = 7!/(7-s)! * C(7-s, m-s) * C(7-m, m-s),
 * so that  d(H) - sum_b <Q_b, M_b(H)> = (144 e(H) M^2 - Zsum(H)) / (5040 M^2)   with Q = Qint / M^2.
 * No canonical root labelling, no automorphism tables: every ordered theta is enumerated.
 *
 * Modes:
 *   rv_eval BLOB list MASKS.bin    MASKS.bin = uint64 producer masks (colex order), e.g. the .npy payload
 *   rv_eval BLOB raw  REPS.bin     REPS.bin = int64 [nreps] own-encoded 6-vertex graphs (4-sets inside {0..5});
 *                                  scans ALL 2^20 links of vertex 6 for every rep (completeness by heredity).
 * Output (stdout): totals, exact minimum numerator Nmin = min(144 e M^2 - Zsum), argmin (own + colex masks),
 * per-(e, cdd) group table "G e cdd count maxZsum argmax_colex", and the K lowest-N graphs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

typedef __int128 i128;
typedef long long ll;

#define NV 7
#define NE 35
static int E7[NE][4];
static int colex_of_own[NE];   /* own index -> colex index */
static int own_of_colex[NE];

typedef struct {
    ll s, m, nsig, nE, conf, nkeys;
    ll *V;            /* [2^nsig] -> key id or -1 */
    ll *T;            /* [2^nE] -> flag index or -1 */
    int nth, nU, npairs;
    int *rootpos;     /* [nth][nsig] */
    int *pos;         /* [nth][nU][nE] (bits nsig.. used) */
    int pairs[64][2];
    int active;
} Block;

static int nblocks, nkeys;
static ll P, Mfac;
static Block B[8];
static ll *Qk[64];
static ll Qn[64];
static int Qnz[64];
static int pmask_n;
static uint64_t pmask[64];
static int cddpairs[70][2];

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
}

static uint64_t colex_to_own(uint64_t h) { uint64_t o = 0; for (int i = 0; i < NE; i++) if ((h >> i) & 1) o |= 1ULL << own_of_colex[i]; return o; }
static uint64_t own_to_colex(uint64_t h) { uint64_t o = 0; for (int i = 0; i < NE; i++) if ((h >> i) & 1) o |= 1ULL << colex_of_own[i]; return o; }

static void rd(void *p, size_t sz, size_t n, FILE *f) { if (fread(p, sz, n, f) != n) { fprintf(stderr, "short read\n"); exit(1); } }

static void load_blob(const char *fn) {
    FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
    ll hdr[4]; rd(hdr, 8, 4, f);
    P = hdr[0]; nblocks = (int)hdr[1]; nkeys = (int)hdr[2]; Mfac = hdr[3];
    for (int b = 0; b < nblocks; b++) {
        Block *X = &B[b];
        ll h[6]; rd(h, 8, 6, f);
        X->s = h[0]; X->m = h[1]; X->nsig = h[2]; X->nE = h[3]; X->conf = h[4]; X->nkeys = h[5];
        ll *order = malloc(sizeof(ll) * X->nE * 4); rd(order, 8, X->nE * 4, f);
        X->V = malloc(sizeof(ll) << X->nsig); rd(X->V, 8, 1ULL << X->nsig, f);
        X->T = malloc(sizeof(ll) << X->nE); rd(X->T, 8, 1ULL << X->nE, f);
        int s = (int)X->s, m = (int)X->m, k = m - s;
        if (2 * m - s > NV || 5040 % X->conf) { fprintf(stderr, "bad block\n"); exit(1); }
        /* ordered root tuples */
        int nth = 1; for (int i = 0; i < s; i++) nth *= (NV - i);
        X->nth = nth;
        X->nU = comb(NV - s, k);
        X->rootpos = malloc(sizeof(int) * nth * (X->nsig ? X->nsig : 1));
        X->pos = malloc(sizeof(int) * nth * X->nU * X->nE);
        int th[8], t = 0;
        /* enumerate ordered s-tuples by odometer over 7^s and keep injective ones */
        int tot = 1; for (int i = 0; i < s; i++) tot *= NV;
        for (int code = 0; code < tot; code++) {
            int c = code, used = 0, ok = 1;
            for (int i = 0; i < s; i++) { th[i] = c % NV; c /= NV; if (used >> th[i] & 1) ok = 0; used |= 1 << th[i]; }
            if (!ok) continue;
            int rest[7], nr = 0; for (int v = 0; v < NV; v++) if (!(used >> v & 1)) rest[nr++] = v;
            for (int j = 0; j < X->nsig; j++)
                X->rootpos[t * X->nsig + j] = own_index(th[order[4*j]], th[order[4*j+1]], th[order[4*j+2]], th[order[4*j+3]]);
            /* U: k-subsets of rest in lex order */
            int u = 0;
            for (int um = 0; um < (1 << nr); um++) {
                if (__builtin_popcount(um) != k) continue;
                int lab[8]; for (int i = 0; i < s; i++) lab[i] = th[i];
                int q = s; for (int i = 0; i < nr; i++) if (um >> i & 1) lab[q++] = rest[i];
                for (int j = 0; j < X->nE; j++)
                    X->pos[(t * X->nU + u) * X->nE + j] = own_index(lab[order[4*j]], lab[order[4*j+1]], lab[order[4*j+2]], lab[order[4*j+3]]);
                u++;
            }
            if (u != X->nU) { fprintf(stderr, "nU mismatch\n"); exit(1); }
            t++;
        }
        if (t != nth) { fprintf(stderr, "nth mismatch\n"); exit(1); }
        /* disjoint ordered pairs of U (by index into the same lex list of subsets of an nr-set) */
        int nr = NV - s, ul[64], nu = 0;
        for (int um = 0; um < (1 << nr); um++) if (__builtin_popcount(um) == k) ul[nu++] = um;
        X->npairs = 0;
        for (int a = 0; a < nu; a++) for (int c = 0; c < nu; c++) if (a != c && !(ul[a] & ul[c])) {
            X->pairs[X->npairs][0] = a; X->pairs[X->npairs][1] = c; X->npairs++;
        }
        /* sanity: conf == nth * npairs */
        if ((ll)nth * X->npairs != X->conf) { fprintf(stderr, "conf mismatch %d %d %lld\n", nth, X->npairs, X->conf); exit(1); }
        free(order);
    }
    for (int i = 0; i < nkeys; i++) {
        ll h[3]; rd(h, 8, 3, f);
        int k = (int)h[0]; Qn[k] = h[1]; Qnz[k] = (int)h[2];
        Qk[k] = malloc(sizeof(ll) * h[1] * h[1]); rd(Qk[k], 8, h[1] * h[1], f);
    }
    fclose(f);
    for (int b = 0; b < nblocks; b++) {
        B[b].active = 0;
        for (ll r = 0; r < (1LL << B[b].nsig); r++) if (B[b].V[r] >= 0 && Qnz[B[b].V[r]]) B[b].active = 1;
    }
    /* p-set masks */
    pmask_n = 0;
    for (int S = 0; S < 128; S++) if (__builtin_popcount(S) == P) {
        uint64_t mk = 0; for (int i = 0; i < NE; i++) { int in = 1; for (int j = 0; j < 4; j++) if (!(S >> E7[i][j] & 1)) in = 0; if (in) mk |= 1ULL << i; }
        pmask[pmask_n++] = mk;
    }
    /* cdd: pairs of edges meeting in exactly one vertex */
    int c = 0;
    for (int i = 0; i < NE; i++) for (int j = i + 1; j < NE; j++) {
        int inter = 0; for (int a = 0; a < 4; a++) for (int b2 = 0; b2 < 4; b2++) inter += E7[i][a] == E7[j][b2];
        if (inter == 1) { cddpairs[c][0] = i; cddpairs[c][1] = j; c++; }
    }
    if (c != 70) { fprintf(stderr, "cdd %d\n", c); exit(1); }
}

static int LAM = 1;   /* admissible: every p-set spans >= LAM edges (env RV_LAM; default 1) */
static inline int admissible(uint64_t h) { for (int i = 0; i < pmask_n; i++) if (__builtin_popcountll(h & pmask[i]) < LAM) return 0; return 1; }

static inline i128 eval_Z(uint64_t h) {
    i128 Zsum = 0;
    int fl[64];
    for (int b = 0; b < nblocks; b++) {
        Block *X = &B[b];
        if (!X->active) continue;
        i128 Zb = 0;
        int nsig = (int)X->nsig, nE = (int)X->nE, nU = X->nU;
        for (int t = 0; t < X->nth; t++) {
            ll r = 0;
            const int *rp = X->rootpos + t * nsig;
            for (int j = 0; j < nsig; j++) r |= (ll)((h >> rp[j]) & 1) << j;
            ll k = X->V[r];
            if (k < 0 || !Qnz[k]) continue;
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
                if (a >= 0 && c >= 0) Zb += Q[a * n + c];
            }
        }
        Zsum += (i128)(5040 / X->conf) * Zb;
    }
    return Zsum;
}

static inline int cdd_of(uint64_t h) { int c = 0; for (int i = 0; i < 70; i++) c += ((h >> cddpairs[i][0]) & (h >> cddpairs[i][1]) & 1); return c; }

static void p128(i128 x, char *buf) {
    char t[64]; int n = 0, neg = x < 0; unsigned __int128 u = neg ? -(unsigned __int128)x : (unsigned __int128)x;
    if (!u) t[n++] = '0';
    while (u) { t[n++] = '0' + (int)(u % 10); u /= 10; }
    int j = 0; if (neg) buf[j++] = '-'; while (n) buf[j++] = t[--n]; buf[j] = 0;
}

#define KLOW 64
typedef struct {
    i128 mz[36][71]; uint64_t ah[36][71]; ll cnt[36][71];
    i128 nmin; uint64_t hmin; ll nadm, nrej;
    i128 lowN[KLOW]; uint64_t lowH[KLOW]; int nlow;
    i128 sumZ; ll sume;
} Acc;

static void acc_init(Acc *A) {
    for (int e = 0; e < 36; e++) for (int c = 0; c < 71; c++) { A->mz[e][c] = 0; A->cnt[e][c] = 0; A->ah[e][c] = 0; }
    A->nmin = 0; A->hmin = 0; A->nadm = 0; A->nrej = 0; A->nlow = 0; A->sumZ = 0; A->sume = 0;
}

static void lowins(Acc *A, i128 N, uint64_t h) {
    if (A->nlow < KLOW) { A->lowN[A->nlow] = N; A->lowH[A->nlow] = h; A->nlow++; }
    else { int w = 0; for (int i = 1; i < KLOW; i++) if (A->lowN[i] > A->lowN[w]) w = i; if (N < A->lowN[w]) { A->lowN[w] = N; A->lowH[w] = h; } }
}

static i128 M2;
static inline void process(Acc *A, uint64_t h) {
    if (!admissible(h)) { A->nrej++; return; }
    A->nadm++;
    i128 Z = eval_Z(h);
    int e = __builtin_popcountll(h), c = cdd_of(h);
    i128 N = (i128)144 * e * M2 - Z;
    if (A->cnt[e][c] == 0 || Z > A->mz[e][c]) { A->mz[e][c] = Z; A->ah[e][c] = h; }
    A->cnt[e][c]++;
    A->sumZ += Z; A->sume += e;
    if (A->nadm == 1 || N < A->nmin) { A->nmin = N; A->hmin = h; }
    lowins(A, N, h);   /* keep the KLOW lowest numerators (low-margin list) */
}

static void merge(Acc *D, Acc *S) {
    for (int e = 0; e < 36; e++) for (int c = 0; c < 71; c++) if (S->cnt[e][c]) {
        if (D->cnt[e][c] == 0 || S->mz[e][c] > D->mz[e][c] || (S->mz[e][c] == D->mz[e][c] && S->ah[e][c] < D->ah[e][c])) { D->mz[e][c] = S->mz[e][c]; D->ah[e][c] = S->ah[e][c]; }
        D->cnt[e][c] += S->cnt[e][c];
    }
    if (S->nadm && (D->nadm == 0 || S->nmin < D->nmin)) { D->nmin = S->nmin; D->hmin = S->hmin; }
    D->nadm += S->nadm; D->nrej += S->nrej; D->sumZ += S->sumZ; D->sume += S->sume;
    for (int i = 0; i < S->nlow; i++) lowins(D, S->lowN[i], S->lowH[i]);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
    init_edges();
    if (getenv("RV_LAM")) LAM = atoi(getenv("RV_LAM"));
    load_blob(argv[1]);
    M2 = (i128)Mfac * Mfac;
    printf("# p=%lld lam=%d nblocks=%d nkeys=%d M=%lld threads=%d\n", P, LAM, nblocks, nkeys, Mfac, omp_get_max_threads());
    for (int b = 0; b < nblocks; b++) printf("# block s=%lld m=%lld conf=%lld nth=%d nU=%d npairs=%d active=%d\n", B[b].s, B[b].m, B[b].conf, B[b].nth, B[b].nU, B[b].npairs, B[b].active);
    fflush(stdout);
    Acc *G = malloc(sizeof(Acc)); acc_init(G);
    double t0 = omp_get_wtime();
    if (!strcmp(argv[2], "list")) {
        FILE *f = fopen(argv[3], "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
        ll n = sz / 8; uint64_t *H = malloc(sz); rd(H, 8, n, f); fclose(f);
        ll off = argc > 4 ? atoll(argv[4]) : 0;   /* header bytes/8 to skip (npy) */
        #pragma omp parallel
        {
            Acc *A = malloc(sizeof(Acc)); acc_init(A);
            #pragma omp for schedule(dynamic, 4096)
            for (ll i = off; i < n; i++) process(A, colex_to_own(H[i]));
            #pragma omp critical
            merge(G, A);
            free(A);
        }
    } else if (!strcmp(argv[2], "raw")) {
        FILE *f = fopen(argv[3], "rb"); fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
        ll nr = sz / 8; uint64_t *Rp = malloc(sz); rd(Rp, 8, nr, f); fclose(f);
        /* link bit j (j-th 3-subset of {0..5}, lex) -> own index of {t, 6} */
        int lpos[20], j = 0;
        for (int a = 0; a < 6; a++) for (int b = a + 1; b < 6; b++) for (int c = b + 1; c < 6; c++) lpos[j++] = own_index(a, b, c, 6);
        uint64_t lo[1024], hi[1024];
        for (int x = 0; x < 1024; x++) { lo[x] = hi[x] = 0; for (int i = 0; i < 10; i++) if (x >> i & 1) { lo[x] |= 1ULL << lpos[i]; hi[x] |= 1ULL << lpos[10 + i]; } }
        for (ll r = 0; r < nr; r++) for (int i = 0; i < NE; i++) if ((Rp[r] >> i & 1) && (E7[i][3] == 6)) { fprintf(stderr, "rep uses vertex 6\n"); return 1; }
        ll ntask = nr * 1024;
        ll done = 0;
        #pragma omp parallel
        {
            Acc *A = malloc(sizeof(Acc)); acc_init(A);
            #pragma omp for schedule(dynamic, 1)
            for (ll task = 0; task < ntask; task++) {
                uint64_t base = Rp[task / 1024] | hi[task % 1024];
                for (int x = 0; x < 1024; x++) process(A, base | lo[x]);
                #pragma omp atomic
                done++;
                if (omp_get_thread_num() == 0 && done % 20480 == 0) { fprintf(stderr, "# progress %lld/%lld %.0fs\n", done, ntask, omp_get_wtime() - t0); }
            }
            #pragma omp critical
            merge(G, A);
            free(A);
        }
    } else { fprintf(stderr, "mode\n"); return 1; }
    double t1 = omp_get_wtime();
    char buf[64], buf2[64];
    p128(G->sumZ, buf); printf("SUMS sumZ %s sume %lld\n", buf, G->sume);
    p128(G->nmin, buf); p128((i128)5040 * M2, buf2);
    printf("TOTAL admissible=%lld rejected=%lld time=%.1fs\n", G->nadm, G->nrej, t1 - t0);
    printf("NMIN %s DEN %s ARGMIN_OWN %llu ARGMIN_COLEX %llu E %d CDD %d\n", buf, buf2, (unsigned long long)G->hmin,
           (unsigned long long)own_to_colex(G->hmin), __builtin_popcountll(G->hmin), cdd_of(G->hmin));
    for (int e = 0; e < 36; e++) for (int c = 0; c < 71; c++) if (G->cnt[e][c]) {
        p128(G->mz[e][c], buf);
        printf("G %d %d %lld %s %llu\n", e, c, G->cnt[e][c], buf, (unsigned long long)own_to_colex(G->ah[e][c]));
    }
    for (int i = 0; i < G->nlow; i++) { p128(G->lowN[i], buf); printf("LOW %s %llu\n", buf, (unsigned long long)own_to_colex(G->lowH[i])); }
    return 0;
}
