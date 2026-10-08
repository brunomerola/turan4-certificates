#!/usr/bin/env bash
# Run the independent checks of this package from its root directory (see README.md and verifier/EXPECTED.txt).
#
#   verifier/verify_all.sh [--fast] [--sigcat] [--lotup] [--raw NAMES] [--producer-n6] [--producer-j4]
#                          [--threads N] [--work DIR]
#
#   --fast          every check except the raw scans (the default when no option is given); about 90 min
#   --sigcat        only the fast checks of parts (m) and (n) (Theorem sigmacat, Remark j4limit); about 6 min
#   --lotup         only part (p) (Proposition lotup: upper bounds for lottery numbers); 1 to 2 min
#   --raw NAMES     the long exhaustive scans, comma-separated: the raw scan of the 4-graph certificates K5_4, K6_4,
#                   K7_4, K5_4minus, K6_4minus (Theorem main), cat_p5_lam3, cat_p6_lam3, cat_p6_lam4, cat_p6_lam5,
#                   cat_p6_lam6, cat_p6_lam9, cat_p6_lam11, cat_p7_lam2, cat_p7_lam3, cat_p7_lam4 (Theorem
#                   catalogue); sigma (Theorem sigma, about 15-20 CPU-min, split over --threads processes); stab_c10
#                   (fact C10 of the stability section over all 604,426 labelled bases, about 20 min); sigcat_J4,
#                   sigcat_K5lt, sigcat_K5_4minus, sigcat_K5eq, sigcat_K5_3minus, sigcat_C5, sigcat_K6_3,
#                   sigcat_K6_4minus, sigcat_K6_4, sigcat_K7_4 (Theorem sigmacat: decoding, then the exact minimum
#                   over all one-vertex extensions; minutes to a few hours each on 2 threads, up to 1.5 GB of memory);
#                   j4_controls (the controls of the dual point of Remark j4limit, 4 to 11 min); or "main" (the
#                   first five), "catalogue" (the ten cat_*), "sigcat" (the ten sigcat_*), "all" (everything).
#                   cat_p5_lam3 and cat_p6_lam11 take seconds, cat_p6_lam9 about 15 min, K5_4minus minutes; each
#                   of the other 4-graph certificate scans 1 to 5 hours on 2 threads (about 5 minutes on 28
#                   threads)
#   --producer-n6   also run the FIRST implementation's checker search/flagalg/verify_cert.py on the three sharp
#                   N = 6 certificates (not an independent check; the independent one is part (j) of --fast)
#   --producer-j4   also run the producer's own checkers search/sigma_catalogue/verify_witness.py and verify_dual3.py
#                   on the two files of Remark j4limit (not independent; the independent checks are part (n))
#   --threads N     OpenMP threads of the C evaluator (default 2)
#   --work DIR      directory for generated files (default: work/ under the package root)
#
# Environment: PYTHON (default python3), CC (default gcc), CFLAGS (default "-O3 -march=native -fopenmp -Wall").
# Every step writes DIR/logs/STEP.log; the summary is DIR/logs/summary.txt. Exit status 0 iff every check passed.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT" || exit 2
PY="${PYTHON:-python3}"
CC="${CC:-gcc}"
CFLAGS="${CFLAGS:--O3 -march=native -fopenmp -Wall}"
FAST=0; SIGCAT=0; LOTUP=0; RAW=""; PRODN6=0; PRODJ4=0; THREADS=2; WORK="work"
while [ $# -gt 0 ]; do
  case "$1" in
    --fast) FAST=1 ;;
    --sigcat) SIGCAT=1 ;;
    --lotup) LOTUP=1 ;;
    --raw) RAW="${2:?--raw needs a list}"; shift ;;
    --producer-n6) PRODN6=1 ;;
    --producer-j4) PRODJ4=1 ;;
    --threads) THREADS="${2:?--threads needs a number}"; shift ;;
    --work) WORK="${2:?--work needs a directory}"; shift ;;
    -h|--help) sed -n '2,31p' "$0"; exit 0 ;;
    *) echo "unknown option: $1 (see --help)"; exit 2 ;;
  esac
  shift
done
if [ "$FAST" = 0 ] && [ "$SIGCAT" = 0 ] && [ "$LOTUP" = 0 ] && [ -z "$RAW" ] && [ "$PRODN6" = 0 ] && [ "$PRODJ4" = 0 ]
then FAST=1; fi
[ "$FAST" = 1 ] && SIGCAT=1      # --fast includes the fast checks of parts (m) and (n)
[ "$FAST" = 1 ] && LOTUP=1       # and part (p)
MAIN5="K5_4 K6_4 K7_4 K5_4minus K6_4minus"
CAT10="cat_p5_lam3 cat_p6_lam3 cat_p6_lam4 cat_p6_lam5 cat_p6_lam6 cat_p6_lam9 cat_p6_lam11 cat_p7_lam2 cat_p7_lam3"
CAT10="$CAT10 cat_p7_lam4"
SIGCAT10="sigcat_J4 sigcat_K5lt sigcat_K5_4minus sigcat_K5eq sigcat_K5_3minus sigcat_C5 sigcat_K6_3"
SIGCAT10="$SIGCAT10 sigcat_K6_4minus sigcat_K6_4 sigcat_K7_4"
RAWLIST=""
for x in ${RAW//,/ }; do
  case "$x" in
    all) RAWLIST="$RAWLIST $MAIN5 $CAT10 sigma stab_c10 $SIGCAT10 j4_controls" ;;
    main) RAWLIST="$RAWLIST $MAIN5" ;;
    catalogue) RAWLIST="$RAWLIST $CAT10" ;;
    sigcat) RAWLIST="$RAWLIST $SIGCAT10" ;;
    *) RAWLIST="$RAWLIST $x" ;;
  esac
done
RAWX=""; RAW4=""      # the special long checks (sigma, stab_c10, sigcat_*, j4_controls) and the 4-graph raw scans
for x in $RAWLIST; do
  case "$x" in
    sigma|stab_c10|j4_controls) RAWX="$RAWX $x" ;;
    sigcat_J4|sigcat_K5lt|sigcat_K5_4minus|sigcat_K5eq|sigcat_K5_3minus|sigcat_C5|sigcat_K6_3) RAWX="$RAWX $x" ;;
    sigcat_K6_4minus|sigcat_K6_4|sigcat_K7_4) RAWX="$RAWX $x" ;;
    *) RAW4="$RAW4 $x" ;;
  esac
done
RAW="${RAWLIST# }"; RAW4="${RAW4# }"; RAWX="${RAWX# }"
case "$WORK" in /*) W="$WORK" ;; *) W="$ROOT/$WORK" ;; esac
V="$ROOT/verifier"
C="$ROOT/certificates"
mkdir -p "$W/logs" || exit 2
SUMMARY="$W/logs/summary.txt"
: > "$SUMMARY"
NPASS=0; NFAIL=0; ELAPSED=0

say() { echo "$*"; echo "$*" >> "$SUMMARY"; }

# run NAME CMD...: run CMD inside the work directory, stdout and stderr to logs/NAME.log; sets ELAPSED, returns rc
run() {
  local name="$1"; shift
  local t0; t0=$(date +%s)
  ( cd "$W" && "$@" ) > "$W/logs/$name.log" 2>&1
  local rc=$?
  ELAPSED=$(( $(date +%s) - t0 ))
  return $rc
}

# runin SUBDIR NAME CMD...: as run, but inside the work subdirectory SUBDIR (created together with SUBDIR/local)
runin() {
  local sub="$1" name="$2"; shift 2
  mkdir -p "$W/$sub/local"
  local t0; t0=$(date +%s)
  ( cd "$W/$sub" && "$@" ) > "$W/logs/$name.log" 2>&1
  local rc=$?
  ELAPSED=$(( $(date +%s) - t0 ))
  return $rc
}

# check NAME RC STRING...: PASS iff RC = 0 and every STRING occurs in logs/NAME.log (fixed-string match);
# a STRING starting with "!" must NOT occur.
check() {
  local name="$1" rc="$2"; shift 2
  local ok=1 why="" s
  [ "$rc" = 0 ] || { ok=0; why="exit status $rc"; }
  for s in "$@"; do
    if [ "${s:0:1}" = "!" ]; then
      grep -qF -- "${s:1}" "$W/logs/$name.log" && { ok=0; why="$why; unexpected: ${s:1}"; }
    else
      grep -qF -- "$s" "$W/logs/$name.log" || { ok=0; why="$why; missing: $s"; }
    fi
  done
  if [ "$ok" = 1 ]; then
    NPASS=$((NPASS + 1)); say "PASS  $name  (${ELAPSED} s)"
  else
    NFAIL=$((NFAIL + 1)); say "FAIL  $name  (${ELAPSED} s)  ${why#; }  -> see $W/logs/$name.log"
  fi
}

# Parameters of the fifteen 4-graph certificates (paper: Theorem main and Table counts; Theorem catalogue and Table
# catalogue): directory, file prefix, p, lambda, 6-vertex representatives, classes on 7 vertices, labelled
# admissible graphs on 7 vertices, admissible raw extensions, exact bound, certificate argmin (colex), raw-scan NMIN
# numerator and denominator 5040 M^2, raw-scan checksum (sumZ, sume; left empty, and then not checked, where the
# independent review did not record it; see EXPECTED.txt).
cinfo() {
  SUMZ=""; SUME=""
  case "$1" in
    K5_4) D=K5_4; PRE=lpcg_p5_conv_full; P=5; LAM=1; NREP=122; NCLS=3908438; NLAB=19199206747; NRAW=86952880
          FRAC=35604499940047/115448720916480; ARG=27303369217; NMIN=1709015997122256; DEN=5541538603991040
          SUMZ=3261680919295407098048; SUME=1670012489 ;;
    K6_4) D=K6_4; PRE=c6l1; P=6; LAM=1; NREP=155; NCLS=7011184; NLAB=34352419335; NRAW=162481445
          FRAC=9224492822869/65970697666560; ARG=2290237506; NMIN=3099429588483984; DEN=22166154415964160
          SUMZ=48794814971023106859848; SUME=2851672393 ;;
    K7_4) D=K7_4; PRE=lpcg_p7_full; P=7; LAM=1; NREP=156; NCLS=7013319; NLAB=34359738367; NRAW=163577855
          FRAC=7476698908057/115448720916480; ARG=33940849153; NMIN=358881547586736; DEN=5541538603991040
          SUMZ=97827306173207768161920; SUME=2862612480 ;;
    K5_4minus) D=K5_4minus; PRE=c5l2_f; P=5; LAM=2; NREP=62; NCLS=360127; NLAB=1741143596; NRAW=10058620
          FRAC=48641756816487/87960930222080; ARG=30440134655; NMIN=12257722717754724; DEN=22166154415964160
          SUMZ=6539164634916749040654; SUME=230039592 ;;
    K6_4minus) D=K6_4minus; PRE=c6l2_f; P=6; LAM=2; NREP=154; NCLS=6986573; NLAB=34244802014; NRAW=160859017
          FRAC=202684381843379/923589767331840; ARG=17444171808; NMIN=4864425164241096; DEN=22166154415964160
          SUMZ=512099105996968813752466; SUME=2833905684 ;;
    cat_p5_lam3) D=catalogue/p5_lam3; PRE=c5l3; P=5; LAM=3; NREP=19; NCLS=3065; NLAB=12420940; NRAW=172932
          FRAC=3298534811415/4398046511104; ARG=34359738367; NMIN=66498461798126400; DEN=88664617663856640
          SUMZ=-291501819268512407200; SUME=4750792 ;;
    cat_p6_lam3) D=catalogue/p6_lam3; PRE=c6l3_f; P=6; LAM=3; NREP=152; NCLS=6836416; NLAB=33540083069
          NRAW=155614014; FRAC=11166973764957/38482906972160; ARG=19736588288; NMIN=6432176888615232
          DEN=22166154415964160 ;;
    cat_p6_lam4) D=catalogue/p6_lam4; PRE=c6l4_f; P=6; LAM=4; NREP=147; NCLS=6292897; NLAB=30911678040
          NRAW=140636249; FRAC=391547931101/1099511627776; ARG=34359738367; NMIN=7893606290996160
          DEN=22166154415964160 ;;
    cat_p6_lam5) D=catalogue/p6_lam5; PRE=c6l5_f; P=6; LAM=5; NREP=138; NCLS=5059432; NLAB=24868266169
          NRAW=112218914; FRAC=64873027783677/153931627888640; ARG=31803041519; NMIN=9341716000849488
          DEN=22166154415964160 ;;
    cat_p6_lam6) D=catalogue/p6_lam6; PRE=c6l6_f; P=6; LAM=6; NREP=123; NCLS=3263333; NLAB=16021163016
          NRAW=74015441; FRAC=612132591153137/1231453023109120; ARG=34342977567; NMIN=44073546563025864
          DEN=88664617663856640; SUMZ=80465255696532640169764; SUME=1490258748 ;;
    cat_p6_lam9) D=catalogue/p6_lam9; PRE=c6l9_f; P=6; LAM=9; NREP=54; NCLS=122149; NLAB=574200468; NRAW=4517328
          FRAC=217134104143267/307863255777280; ARG=4294711295; NMIN=15633655498315224; DEN=22166154415964160
          SUMZ=-4041155699855865843576; SUME=112293269 ;;
    cat_p6_lam11) D=catalogue/p6_lam11; PRE=c6l11_f; P=6; LAM=11; NREP=18; NCLS=1939; NLAB=7009743; NRAW=155364
          FRAC=24267758577011/28862180229120; ARG=2146992127; NMIN=1164852411696528; DEN=1385384650997760
          SUMZ=-16854299946483882996; SUME=4433920 ;;
    cat_p7_lam2) D=catalogue/p7_lam2; PRE=c7l2_f; P=7; LAM=2; NREP=156; NCLS=7013318; NLAB=34359738332
          NRAW=163577834; FRAC=111867473425/1099511627776; ARG=34359738367; NMIN=563812066062000
          DEN=5541538603991040 ;;
    cat_p7_lam3) D=catalogue/p7_lam3; PRE=c7l3_f; P=7; LAM=3; NREP=156; NCLS=7013315; NLAB=34359737737
          NRAW=163577622; FRAC=4043518541347/30786325577728; ARG=30064771072; NMIN=727833337442460
          DEN=5541538603991040 ;;
    cat_p7_lam4) D=catalogue/p7_lam4; PRE=c7l4_f; P=7; LAM=4; NREP=156; NCLS=7013305; NLAB=34359731192
          NRAW=163576247; FRAC=11199502214035/69269232549888; ARG=7918845952; NMIN=3583840708491200
          DEN=22166154415964160 ;;
    *) echo "unknown certificate name: $1 (see --help)"; exit 2 ;;
  esac
  if [ "$LAM" = 1 ]; then SUF=""; else SUF="_l$LAM"; fi
}
for n in $RAW4; do cinfo "$n"; done      # reject unknown names before any work

say "# verify_all.sh  $(date -u +%FT%TZ)  root=$ROOT  work=$W  threads=$THREADS  python=$("$PY" -c 'import sys; print(sys.version.split()[0])' 2>&1)"
say "# options: fast=$FAST sigcat=$SIGCAT lotup=$LOTUP raw='${RAW}' producer_n6=$PRODN6 producer_j4=$PRODJ4"

# --- 0. integrity and build -----------------------------------------------------------------------------------------
if command -v sha256sum > /dev/null; then SHA="sha256sum"; else SHA="shasum -a 256"; fi
t0=$(date +%s); (cd "$ROOT" && $SHA -c SHA256SUMS) > "$W/logs/integrity.log" 2>&1; rc=$?; ELAPSED=$(( $(date +%s) - t0 ))
check integrity $rc "!FAILED"
run npz_hashes "$PY" -c '
import hashlib, json, sys, glob, os
ok = True
js = sorted(glob.glob(os.path.join(sys.argv[1], "**", "*.cert.json"), recursive=True))
for j in js:
    d = json.load(open(j))
    h = hashlib.sha256(open(j[:-5] + ".npz", "rb").read()).hexdigest()
    print(os.path.relpath(j, sys.argv[1]), "cert_npz_sha256", d["cert_npz_sha256"] == h)
    ok &= d["cert_npz_sha256"] == h
print(len(js), "certificate files")
print("ALL NPZ HASHES MATCH" if ok else "NPZ HASH MISMATCH")' "$C"
check npz_hashes $? "26 certificate files" "ALL NPZ HASHES MATCH"
NEED_EVAL=0
{ [ "$FAST" = 1 ] || [ -n "$RAW4" ]; } && NEED_EVAL=1
if [ "$NEED_EVAL" = 1 ]; then
  # shellcheck disable=SC2086
  run build "$CC" $CFLAGS -o rv_eval "$V/rv_eval.c"
  check build $?
fi

# --- (a) and (h): the 4-graph certificates (five of Theorem main, ten of Theorem catalogue) -------------------------
PREPARED=" "
prepare() {   # own 6-vertex representatives and Burnside counts, then the evaluation tables of the certificate
  local n="$1"
  case "$PREPARED" in *" $n "*) return ;; esac
  PREPARED="$PREPARED$n "
  cinfo "$n"
  run "reps_$n" "$PY" "$V/rv_reps.py" "$P" "$LAM"
  check "reps_$n" $? '"burnside_classes"'
  run "counts_$n" "$PY" -c '
import json, sys
d = json.load(open(sys.argv[1]))
got = [d["n_reps6"], d["burnside_classes"], d["labelled_7"], d["raw_extensions"]]
exp = [int(x) for x in sys.argv[2:]]
print("6-vertex reps, classes on 7 vertices (Burnside), labelled admissible, raw extensions:", got, "expected", exp)
print("COUNTS OK" if got == exp else "COUNTS DIFFER")' "reps_p${P}${SUF}.json" "$NREP" "$NCLS" "$NLAB" "$NRAW"
  check "counts_$n" $? "COUNTS OK"
  run "prep_$n" env RV_LAM="$LAM" "$PY" "$V/rv_prep.py" "$C/$D/$PRE" "blob_$n.bin"
  check "prep_$n" $?
  run "lists_$n" "$PY" -c '
import json, sys
r = json.load(open(sys.argv[1]))
ks = [k for b in r["blocks"] for k in b["keys"]]
bad = [k["k"] for k in ks if k["missing"] or k["inadmissible_listed"] or k["dup_classes"]]
print("keys", len(ks), "; factor rows per key", [k["rows"] for k in ks])
print("flag lists: root part = sigma (asserted by rv_prep), pairwise non-isomorphic, complete, all admissible:",
      "LISTS OK" if not bad else "LISTS BAD %s" % bad)' "blob_$n.bin.json"
  check "lists_$n" $? "LISTS OK"
}

if [ "$FAST" = 1 ]; then
  for n in $MAIN5 $CAT10; do
    prepare "$n"; cinfo "$n"
    run "pycheck_cert_$n" env RV_LAM="$LAM" "$PY" "$V/rv_pycheck.py" "$C/$D/$PRE" "$ARG"
    check "pycheck_cert_$n" $? "$ARG $FRAC " "== cert bound: True"
  done
  # (h) the exact bounds of the catalogue and their decimals in Table catalogue of the paper (t floored, pi ceiled)
  run cat_decimals "$PY" -c '
import json, os, sys
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
getcontext().prec = 80
def dec(x, mode):
    return str((Decimal(x.numerator) / Decimal(x.denominator)).quantize(Decimal("1e-12"), rounding=mode))
ok = True
for row in sys.argv[2:]:
    d, pre, frac, t12, pi12 = row.split(":")
    c = json.load(open(os.path.join(sys.argv[1], d, pre + ".cert.json")))
    b, p = Fraction(c["bound"]), c["p"]
    tl, pu = dec(b, ROUND_FLOOR), dec(1 - b, ROUND_CEILING)
    good = b == Fraction(frac) and tl == t12 and pu == pi12
    ok &= good
    print(d, "p", p, "bound", b, "as expected:", b == Fraction(frac), "; t >", tl, "(table:", t12 + ") ; pi <", pu,
          "(expected:", pi12 + ") ;", good)
print("CATALOGUE DECIMALS OK" if ok else "CATALOGUE DECIMALS DIFFER")' "$C/catalogue" \
    p5_lam3:c5l3:3298534811415/4398046511104:0.749999983648:0.250000016352 \
    p6_lam3:c6l3_f:11166973764957/38482906972160:0.290180099259:0.709819900741 \
    p6_lam4:c6l4_f:391547931101/1099511627776:0.356110768826:0.643889231174 \
    p6_lam5:c6l5_f:64873027783677/153931627888640:0.421440536123:0.578559463877 \
    p6_lam6:c6l6_f:612132591153137/1231453023109120:0.497081561103:0.502918438897 \
    p6_lam9:c6l9_f:217134104143267/307863255777280:0.705293990330:0.294706009670 \
    p6_lam11:c6l11_f:24267758577011/28862180229120:0.840815156178:0.159184843822 \
    p7_lam2:c7l2_f:111867473425/1099511627776:0.101742874380:0.898257125620 \
    p7_lam3:c7l3_f:4043518541347/30786325577728:0.131341381781:0.868658618219 \
    p7_lam4:c7l4_f:11199502214035/69269232549888:0.161680760732:0.838319239268
  check cat_decimals $? "CATALOGUE DECIMALS OK"
fi
for n in $RAW4; do
  prepare "$n"; cinfo "$n"
  run "raw_$n" env RV_LAM="$LAM" OMP_NUM_THREADS="$THREADS" ./rv_eval "blob_$n.bin" raw "reps_p${P}${SUF}.bin"
  rc=$?
  if [ -n "$SUMZ" ]; then
    check "raw_$n" $rc "TOTAL admissible=$NRAW " "NMIN $NMIN DEN $DEN " "SUMS sumZ $SUMZ sume $SUME"
  else
    check "raw_$n" $rc "TOTAL admissible=$NRAW " "NMIN $NMIN DEN $DEN "
  fi
  run "compare_raw_$n" env RV_LAM="$LAM" "$PY" "$V/rv_compare.py" "$C/$D/$PRE.cert.json" "logs/raw_$n.log" --raw
  check "compare_raw_$n" $? "rv bound $FRAC = " "EQUAL: True" "same set: True; maxZ equal on all: True"
  RARG=$(awk '$1 == "NMIN" {print $8}' "$W/logs/raw_$n.log")
  run "pycheck_raw_$n" env RV_LAM="$LAM" "$PY" "$V/rv_pycheck.py" "$C/$D/$PRE" "${RARG:-0}"
  check "pycheck_raw_$n" $? "${RARG:-0} $FRAC " "== cert bound: True"
done

# --- (k) and (l): sigma(K_5^(4)) and the stability facts; helpers used by --fast and by --raw sigma / stab_c10 -------
SIGMA_MASKS="2290237506 17418883036 19736566032 26383237633 26894418433 29337331976 31813832752 31837153552"
SIGMA_MASKS="$SIGMA_MASKS 34326183967 34359738367 26887858440 31836891137 31830626344 22674440848 239042524 31826919937"
SIGMA_MASKS="$SIGMA_MASKS 31446477202 23327952978 4413698879 18972461960 33940307901 19616300926 25201702534 17179831775"
SIGMA_MASKS="$SIGMA_MASKS 32305321325 23724464310"
MU0=176165518826891630107/6280747422216628134215680
SIGPREP=0
sigma_prepare() {   # build the two evaluators, own 6-vertex representatives, exact Gram matrices and tables (once)
  [ "$SIGPREP" = 1 ] && return
  SIGPREP=1
  mkdir -p "$W/sigma/local"
  # shellcheck disable=SC2086
  run sigma_build "$CC" $CFLAGS -DNLC=4 -DTWOPHASE -o sigma/rs_eval4 "$V/rs_eval.c"
  rc=$?
  # shellcheck disable=SC2086
  [ "$rc" = 0 ] && { ( cd "$W" && "$CC" $CFLAGS -DNLC=4 -DSYM -o sigma/rs_eval4s "$V/rs_eval.c" ) >> "$W/logs/sigma_build.log" 2>&1; rc=$?; }
  check sigma_build $rc
  runin sigma sigma_reps "$PY" "$V/rs_reps.py"
  check sigma_reps $? '"admissible_6vertex_labelled": 27449' '"reps6": 122' '"raw_admissible_extensions": 86952880' \
    '"labelled_admissible_7vertex": 19199206747'
  runin sigma sigma_prep "$PY" "$V/rs_prep.py" --sharp "$C/sigma_K5_4/sharp_klp5.cert.npz" \
    --keys "$C/sigma_K5_4/sharp_klp5.cert.json" --keys-check "$C/K5_4/lpcg_p5_conv_full.cert.json" --target 33/64 \
    --out local/blob_sharp.bin --sym-out local/blob_sharp_sym.bin
  check sigma_prep $? '"k": 6, "W_rows": 272, "B_rows": 20, "B_max": 16, "X_dim": 20, "X_psd_exact": true' \
    '"k": 8, "W_rows": 160, "B_rows": 5, "B_max": 2, "X_dim": 5, "X_psd_exact": true' \
    '"k": 10, "W_rows": 106, "B_rows": 1, "B_max": 1, "X_dim": 1, "X_psd_exact": true' \
    '"n_flags": 1024, "distinct_classes": 1024, "all_admissible": true, "admissible_classes": 1024, "missing": 0, "extra": 0' \
    '"n_flags": 809, "distinct_classes": 809, "all_admissible": true, "admissible_classes": 809, "missing": 0, "extra": 0' \
    "D 325412473936958852881289201034484891095266230272 limbs 4 G_max_bits 163 unused [] sym" \
    "'aut_sizes': {0: 1, 1: 2, 2: 6, 3: 6, 4: 24, 5: 24, 6: 24, 7: 12, 8: 12, 9: 24, 10: 120}" \
    '!"missing": 1' '!"all_admissible": false'
}
STABPREP=0
stab_prepare() {    # labelled Giraud systems on 7, 8, 9 vertices and their classes (once)
  [ "$STABPREP" = 1 ] && return
  STABPREP=1
  runin stab stab_cf "$PY" "$V/r_cf.py"
  check stab_cf $? "'f5_equals_definition': True, 'giraud_tests_agree': True" \
    '"labelled_normalised": 2404, "labelled_all_matrices": 2404, "same_set": true, "partition_sets_equal": true, "classes": 10, "orbit_sizes": [1, 7, 21, 30, 35, 105, 210, 315, 420, 1260], "all_odd": true, "tight_masks_are_giraud": true, "each_tight_mask_in_exactly_one_orbit_and_bijective": true, "multi_partition": {"7": 30}' \
    '"part_sizes": [[3, 4], [3, 4], [3, 4], [3, 4], [3, 4], [3, 4], [3, 4]], "five_counts_hist": [0, 21, 0, 0, 0, 0], "complements_of_edges_form_Fano_plane": true, "orbit_size": 30, "orbit_contains_tight_mask": [2290237506]' \
    '"labelled_normalised": 32981, "labelled_all_matrices": 32981, "same_set": true, "partition_sets_equal": true, "classes": 22,' \
    '"multi_partition": {"7": 30}, "multi_partition_orbits": [{"rep": 42667670269382951568, "orbit_size": 30, "edges": 14, "partitions": 7,' \
    '"partition_parts_are_blocks": true, "five_counts_hist": [0, 56, 0, 0, 0, 0], "steiner_quadruple_system_S(3,4,8)": true, "blocks_closed_under_complement": true, "all_multi_partition_graphs_in_this_orbit": true' \
    '"labelled_normalised": 604426, "classes": 33,' '"orbit_sizes_sum": 604426, "multi_partition": {}, "odd_on_2000_sample": true'
}

if [ "$FAST" = 1 ]; then
  # --- (b) the three 5-graph certificates -------------------------------------------------------------------------
  for x in "k65 888 5410603883828909241569858123831800227/23261489926236027775816623554906030080" \
           "h45 150 28392597905796184739036021237749789573/46522979852472055551633247109812060160" \
           "k75 1043 1006048490611906127582298374928111023/13956893955741616665489974132943618048"; do
    set -- $x
    run "r5_$1" "$PY" "$V/rv_r5.py" "$C/fivegraphs/cert_$1_N7_all.json"
    check "r5_$1" $? "exact min over all $2 classes: $3 = " "EQUAL: True; PSD: True"
  done

  # --- (c) the N = 6 dual points (Theorem n6) ---------------------------------------------------------------------
  for x in "p5 5 1 1/4" "p6 6 1 1/10" "p7 7 1 0" "p5_l2 5 2 1/2"; do
    set -- $x
    run "dual6_$1" "$PY" "$V/rv_dual6.py" "$C/n6_dual_points/n6_control_$1.json" "$2" "$3"
    check "dual6_$1" $? "V = sum y d = $4; claimed $4; equal True" "10 moment matrices; ALL PSD: True"
  done

  # --- (d) the N = 7 dual point (Theorem n7opt) -------------------------------------------------------------------
  run selftest "$PY" "$V/rv1_selftest.py"
  check selftest $? "K_7 sanity OK" "!False"
  run dual7 "$PY" "$V/rv1_dual.py" "$C/n7_dual_point/dual_Wodd_best.json" "$C/K5_4/lpcg_p5_conv_full" rv1_dual.json
  check dual7 $? "sum num == den: True all num > 0: True" \
    "pairwise non-isomorphic (canonical forms over all 5040 permutations): True" \
    "admissible (every 5-set has an edge): True  odd (every 5-set spans 1,3,5 edges): True" \
    "V == JSON value: True  V == 1336237682928914994292138923/(7*2^89): True" \
    "ALL moment matrices PSD (all labelled types, all six blocks): True" \
    "negative control all PSD: False (expected False)" \
    "b_cert <= sum y (d - c) <= V : True" \
    "claimed 0.308401201195 >= V: True  claimed 0.691598798805 <= 1-V: True"
  run lift "$PY" "$V/rv1_lift.py" "$C/n7_dual_point/dual_Wodd_best.json" 30
  check lift $? "identities checked 17045; all equal: True"

  # --- (e) Giraud's construction against the t(5,4) certificate, and the parity remark ----------------------------
  run giraud "$PY" "$V/rv1_giraud.py" "$C/K5_4/lpcg_p5_conv_full" giraud
  check giraud $? "total probability == 1: True ; distinct labelled graphs: 2404" "all odd: True ; all admissible: True" \
    "E[d] = 5/16 == 5/16: True" "isomorphism classes: 10 " \
    "value constant on every isomorphism class (evaluator self-consistency): True" \
    "claimed 1091855/15393162788864: True" "claimed 52186103/115448720916480: True" \
    "claimed 173614396936021/2^49: True" "claimed 2307463508139/2^49: True" "all support graphs >= b: True" \
    "Giraud law: all moment matrices PSD: True"
  run claim3_AB "$PY" "$V/rv1_claim3.py" "$C/K5_4/lpcg_p5_conv_full" claim3 AB
  check claim3_AB $? "even 5-sets: 14 " "own evaluator: d - c = 35604499940047/115448720916480  == b: True" \
    "generated: 1048576  distinct: True" "all odd: True" "dual-point support graphs contained: True" \
    "Giraud-law graphs contained: True" "isomorphism classes by Burnside: 423"
  run odd_scan env OMP_NUM_THREADS="$THREADS" ./rv_eval blob_K5_4.bin list claim3/odd_all.u64
  check odd_scan $? "TOTAL admissible=1048576 rejected=0 "
  run claim3_C "$PY" "$V/rv1_claim3.py" "$C/K5_4/lpcg_p5_conv_full" claim3 C logs/odd_scan.log
  check claim3_C $? "claimed 376193/28862180229120: True" "own == C numerator: True"

  # --- (f) r = 3 comparison: plain flag algebras on 5 vertices for K_4^(3) ----------------------------------------
  run k43_n5 "$PY" "$V/indep_k43_n5.py"
  check k43_n5 $? "labelled admissible graphs 768" "exact b > Chung-Lu: True"
  if grep -qF "exact b = 60084175574225/140737488355328 " "$W/logs/k43_n5.log"; then
    say "INFO  k43_n5: exact b = 60084175574225/140737488355328 (the value quoted in the paper)"
  else
    say "INFO  k43_n5: exact b differs from the paper's 60084175574225/2^47 (it depends on the floating-point SDP" \
        "solution, i.e. on the clarabel version); the check above only requires b > Chung-Lu. Got: $(grep 'exact b =' "$W/logs/k43_n5.log")"
  fi

  # --- (i) the seven-vertex limits of the other certificates (Theorem n7opt, Table n7limits) ----------------------
  mkdir -p "$W/results"
  run dual2_selftest "$PY" "$V/rv2_selftest.py"
  check dual2_selftest $? \
    "flag_universe_pl(p=5, lam=1) == rv1_lib.flag_universe for every labelled type of all six blocks" \
    "canonical_forms == rv1_lib.canonical_graph on 80 graphs (incl. relabelled copies); admissible() == brute force" \
    "self-tests OK"
  # tag | support size | V | (5,6) summary: (occurring, live, PSD) by e(type) | ceil12(V) | floor12(1 - V)
  for x in "p6_f3|884|37867960095539974857512229/270799383593676935134183424|{0: [(371, 371, True)], 1: [(221, 221, True)], 2: [(74, 74, True)], 3: [(42, 42, True)], 4: [(40, 40, True)], 5: [(237, 237, True)]}|0.139837689411|0.860162310589" \
           "p7_f1|1058|6413961863433425283788748011/99035203142830421991929937920|{0: [(422, 422, True)], 1: [(269, 269, True)], 2: [(104, 104, True)], 3: [(63, 63, True)], 4: [(97, 97, True)], 5: [(217, 217, True)]}|0.064764464150|0.935235535850" \
           "h44_nt3|1583|384021679403583125317285924719/693246421999812953943509565440|{2: [(373, 373, True)], 3: [(457, 457, True)], 4: [(588, 588, True)], 5: [(768, 768, True)]}|0.553946861054|0.446053138946" \
           "k6m_W|2126|153033054540662097928475789999/693246421999812953943509565440|{0: [(371, 371, True)], 1: [(417, 417, True)], 2: [(233, 233, True)], 3: [(207, 207, True)], 4: [(187, 187, True)], 5: [(454, 454, True)]}|0.220748423194|0.779251576806" \
           "p5no56_W|1569|390498160148386371476706381487/1386492843999625907887019130880||0.281644555065|0.718355444935" \
           "p5no56_Wodd|55|24406135009720244409729340413/86655802749976619242938695680||0.281644555070|0.718355444930"; do
    IFS='|' read -r tag nsup val s56 c12 f12 <<< "$x"
    run "dual2_$tag" "$PY" "$V/rv2_dual.py" "$tag" "$C/n7_dual_points/dual_$tag.json"
    rc=$?
    common=("JSON p/lam/blocks equal the target: True"
            "{'support': $nsup, 'den_is_2^110': True, 'sum_is_den': True, 'all_num_ge_1': True, 'masks_in_range': True, 'distinct_masks': True, 'pairwise_non_isomorphic': True, 'admissible': True}"
            "{'V': '$val', " "'V_eq_json': True, 'V_eq_claim': True}"
            "ALL listed blocks PSD (every labelled type): True"
            "A-2lam xx^T exactly indefinite=True, residual PASS=False (expected False); A-lam/2 xx^T residual PASS=True (expected True)"
            "\"ceil12_V\": \"$c12\", \"floor12_1-V\": \"$f12\", \"ceil12_claim_ok\": true, \"pifloor_claim_ok\": true"
            "OVERALL $tag: PASS")
    case "$tag" in
      p5no56_W) check "dual2_$tag" $rc "${common[@]}" \
          "CONTROL (5,6) for the no-(5,6) point: 31 of 31 occurring labelled types FAIL (expected >= 1)" \
          "\"b<=sum\": true" "\"min_support_ge_b\": true" ;;
      p5no56_Wodd) check "dual2_$tag" $rc "${common[@]}" \
          "CONTROL (5,6) for the no-(5,6) point: 16 of 16 occurring labelled types FAIL (expected >= 1)" ;;
      *) check "dual2_$tag" $rc "${common[@]}" \
          "summary (5, 6): e(type) -> sorted set of (occurring, live, PSD): $s56" \
          "\"b<=sum\": true, \"sum<=V\": true" "\"min_support_ge_b\": true" "\"all_key_QY_ge_0\": true" \
          "\"b<=V\": true" "\"picert_claim_ok\": true" "\"b_eq_cert_bound\": true" ;;
    esac
  done
  run dual2_lift "$PY" "$V/rv2_lift.py" 40 20
  check dual2_lift $? "identities checked 27254; all equal: True"
  run dual2_decimals "$PY" "$V/rv2_decimals.py"
  check dual2_decimals $? "ALL OK" "!BAD"

  # --- (j) the sharp N = 6 certificates (Remark n6sharp), independent checker -------------------------------------
  for x in "p5 1/4 122 27449" "p6 1/10 155 32767" "p5_l2 1/2 62 12068"; do
    set -- $x
    run "sharp6_$1" "$PY" "$V/rv_sharp6.py" "$C/n6_dual_points/cert_turan_r4_$1_N6_all_sharp.json"
    check "sharp6_$1" $? "all symmetric: True; all PSD (exact L D L^T, product checked): True" \
      "flags: all sigma-flags, pairwise non-isomorphic: True" \
      "labelled admissible 4-graphs on 6 vertices: $4; isomorphism classes: $3 (orbit sizes sum to the labelled count: True)" \
      "exact min over all $3 classes: b = $2 = " "random relabelling of every class equals its value): True" \
      "claimed $2; EQUAL: True" "SHARP N = 6 CERTIFICATE OK"
  done
  run sharp6_negative "$PY" -c '
import json, subprocess, sys
chk, src = sys.argv[1], json.load(open(sys.argv[2]))
def variant(name, f):
    d = json.loads(json.dumps(src))
    f(d)
    fn = "sharp6_negative_" + name + ".json"
    json.dump(d, open(fn, "w"))
    r = subprocess.run([sys.executable, chk, fn], capture_output=True, text=True)
    rej = r.returncode == 1 and "SHARP N = 6 CERTIFICATE NOT OK" in r.stdout
    print(name, "rejected:", rej)
    return rej
def b45(d):
    return [b for b in d["blocks"] if b["s"] == 4][0]
def claim(d):
    d["bound"] = "1/3"
def notpsd(d):
    b45(d)["Qnum"][0][0] = "-1"
def double(d):
    b45(d)["den"] = str(int(b45(d)["den"]) // 2)
ok = all([variant("claimed_bound_1_3", claim), variant("Q_not_psd", notpsd), variant("Q_doubled", double)])
print("NEGATIVE CONTROLS REJECTED" if ok else "A NEGATIVE CONTROL WAS ACCEPTED")' \
    "$V/rv_sharp6.py" "$C/n6_dual_points/cert_turan_r4_p5_N6_all_sharp.json"
  check sharp6_negative $? "NEGATIVE CONTROLS REJECTED"
  # --- (k) sigma(K_5^(4)) = 31/64 (Theorem sigma): the sharp certificate --------------------------------------------
  run sigma_keys "$PY" -c '
import hashlib, json, sys
import numpy as np
a, b, npz = json.load(open(sys.argv[1])), json.load(open(sys.argv[2])), sys.argv[3]
z = np.load(npz)
same = a["keys"] == b["keys"]
nums = int(z["kb"]) == a["kb"] and str(z["q"]) == a["q"] and str(z["D"]) == a["D"] and int(a["D"]) == 57 * 2 ** 152
sha = hashlib.sha256(open(npz, "rb").read()).hexdigest() == a["cert_npz_sha256"]
shape = all(z["W%d" % k].shape[1] == key["n_flags"] == len(key["flags"]) for k, key in enumerate(a["keys"]))
print("keys", len(a["keys"]), "identical to the K5_4 key list:", same, "; kb, q, D = 57*2^152 as in the npz:", nums,
      "; npz sha256:", sha, "; W widths = flag counts:", shape)
print("SIGMA KEY LIST OK" if same and nums and sha and shape else "SIGMA KEY LIST BAD")' \
    "$C/sigma_K5_4/sharp_klp5.cert.json" "$C/K5_4/lpcg_p5_conv_full.cert.json" "$C/sigma_K5_4/sharp_klp5.cert.npz"
  check sigma_keys $? "SIGMA KEY LIST OK"
  runin sigma sigma_reduction "$PY" "$V/rs_reduction.py"
  check sigma_reduction $? "(a) co2 / s_n / f identities on 52 random graphs: True" \
    "(b) E_S f(H[S]) == f_n(H) exactly on 9 graphs (n = 8, 9, 10): True" "(c) K5-free(G) <=> admissible(H): True" \
    "(d) Giraud exact limit: s(G) = 31/64  d(G) = 11/16" "ALL OK"
  runin sigma sigma_giraud7 "$PY" "$V/rs_giraud7.py"
  check sigma_giraud7 $? '"labelled_graphs_in_support": 2404' '"outside_the_ten_orbits": 0' '"E_phi_f": "33/64"' '"PASS": true'
  sigma_prepare
  # shellcheck disable=SC2086
  runin sigma sigma_pycheck "$PY" "$V/rs_pycheck.py" "$C/sigma_K5_4/sharp_klp5.cert.npz" \
    "$C/sigma_K5_4/sharp_klp5.cert.json" 33/64 $SIGMA_MASKS
  set --
  for x in 2290237506 17418883036 19736566032 26383237633 26894418433 29337331976 31813832752 31837153552 34326183967 \
           34359738367; do set -- "$@" "PY $x slack 33/64 slack-b 0 S 0"; done
  check sigma_pycheck $? "$@" "PY 26887858440 slack 3238686555099275773335067/6280747422216628134215680 slack-b $MU0 S 46001745554897451729515273061756040526495219712"
  # shellcheck disable=SC2086
  runin sigma sigma_each_masks "$PY" -c 'import sys, numpy as np; np.array([int(x) for x in sys.argv[1:]], dtype=np.uint64).tofile("local/check_masks.u64"); print(len(sys.argv) - 1, "masks written")' $SIGMA_MASKS
  check sigma_each_masks $? "26 masks written"
  runin sigma sigma_each ./rs_eval4 local/blob_sharp.bin each local/check_masks.u64
  check sigma_each $? "EACH colex 2290237506 F 168 S 0 0 0 0" "EACH colex 34359738367 F 420 S 0 0 0 0" "!inadmissible"
  runin sigma sigma_each_cmp "$PY" "$V/rs_cmp_py.py" "$W/logs/sigma_each.log" "$W/logs/sigma_pycheck.log"
  check sigma_each_cmp $? "26/26 identical"

  # --- (l) stability for sigma(K_5^(4)): the computer facts of the stability section ---------------------------------
  stab_prepare
  runin stab stab_lg8 "$PY" "$V/r_lg.py" 8
  check stab_lg8 $? '"labelled_R": 2404, "survivors_total": 32981, "G_n_size": 32981, "survivors_equal_giraud_extensions_for_every_R": true, "non_giraud_survivors": [], "C8_PASS": true'
  runin stab stab_lg9 "$PY" "$V/r_lg.py" 9
  check stab_lg9 $? '"labelled_R": 32981, "survivors_total": 604426, "G_n_size": 604426, "survivors_equal_giraud_extensions_for_every_R": true, "non_giraud_survivors": [], "C9_PASS": true'
  runin stab stab_lg10 "$PY" "$V/r_lg.py" 10
  check stab_lg10 $? '"locally_giraud_total": 1151, "C10_PASS": true' "!False"
  runin stab stab_margin "$PY" "$V/r_thm1.py"
  check stab_margin $? '"lambda": "1354230560069659129/4722366482869645213696"' '"lambda": "0"' \
    '"lambda": "6614817002453309619559/4722366482869645213696"' '"lambda": "1930498463925983500877/4722366482869645213696"' \
    '"lambda": "1262520681467868487601/2361183241434822606848"' '"lambda": "13141619641562093306051/1180591620717411303424"' \
    '"lambda": "55057280191985960149631/4722366482869645213696"' '"lambda": "48513285016664528863283/2361183241434822606848"' \
    '"lambda": "6429978663646678105565/590295810358705651712"' '"lambda": "16464561510059415945793/2361183241434822606848"' \
    '"p_over_1mp_formulas_ok_10..2000": true' '"analytic_block_bounds_ok_10..20000": true' \
    '"constant_sum_form_le_68": true' '"c(n)<=68/(n-7)_exact_10..20000": true' "\"mu\": \"$MU0\"" \
    '"inv_mu_le_35653": true' '"K_le_7.81e8": true' '"rho_factor": "221/320"' '"120*126/rho_factor_le_21894": true' \
    '"mu/1200_implies_f_le_33/64+1/100": true' '"1200/mu_le_K": true' '"frames_min_over_a": 1200' \
    '!"min_diag_nonneg": false' '!"max_abs_entry_equals_lambda": false'
  runin stab stab_formF "$PY" "$V/r_formF.py"
  check stab_formF $? \
    "7 {'G_{m+1}': 32981, 'distinct_links': 32981, 'generated_formF': 32981, 'generated_equals_links': True, 'lemma1_i_ok': True, 'with_several_reps': 30, 'lemma2_ok': True, 'generated_reps_equal_predicate_reps': True}" \
    "6 {'G_{m+1}': 2404, 'distinct_links': 2404, 'generated_formF': 2404, 'generated_equals_links': True, 'lemma1_i_ok': True, 'with_several_reps': 30, 'lemma2_ok': True, 'generated_reps_equal_predicate_reps': True}" \
    "5 {'G_{m+1}': 197, 'distinct_links': 197, 'generated_formF': 197, 'generated_equals_links': True, 'lemma1_i_ok': True, 'with_several_reps': 15, 'lemma2_ok': True, 'generated_reps_equal_predicate_reps': True, 'all_3graphs_formF_equals_generated': True}" \
    "!False"
  runin stab stab_cor "$PY" "$V/r_cor.py"
  check stab_cor $? "{'instances': 174, 'formula_equals_bruteforce': True}" "'expansion_identities_ok': True" \
    '"delta_Lambda_removal": "-28/5", "equals_-(n-4)/5": true, "admissible_after": true' \
    "'all_equal': True, 'min_nonzero_over_a': 1200, '(10)_6': 151200" "'violations': 0" "'2d-d^2_at_5/16': '135/256'" \
    '!"cross_removal_characterisation_ok": false' '!"all_other_admissible_flips_ge_(n-4)/10": false' \
    '!"equals_-(n-4)/5": false' '!"admissible_after": false'
fi

# --- (m) and (n): Theorem sigmacat and Remark j4limit (codegree-squared densities, Table sigmacat) ------------------
# The ten certificates of Theorem sigmacat: directory, file prefix (the checkers append .cert.json / .cert.npz), key
# list, problem name of the checkers, the review whose programs check it (s7 = R7_SIG7 for J4, K5lt, K5_4minus; s8 =
# R7_SIG8 for the other seven), prep program and evaluator (the J_4 certificate has M = 2^30 and needs the __int128
# Gram variant), the exact bound b, the stated sigma bound 1 - b, its 12-digit ceiling and the 7-digit ceiling of Table
# sigmacat, keys and flags; from the logs of the independent reviews: 6-vertex representatives, admissible and rejected
# one-vertex extensions, labelled extensions at the minimum, objective groups, labelled admissible 7-vertex graphs, the
# a-priori |Z| bound and the sha256 of the evaluation table written by the prep; the number of spot-check graphs and
# of those at b (sigcat_spot_values.txt).
scinfo() {
  SZB=""; SPC=""; SPROG=s7
  case "$1" in
    J4) SD=sigma_catalogue/J4; SPRE=cosig3_s3_J4_r30.cosig3; SKEYS=cosig3_s3_J4_r30.keys.json; SPROB=s3_J4
        SPREP=s7_prep128.py; SEVAL=s7_eval128; SRN=J4_r30
        SB=23220181141640969669/32281802128991715328; SSIG=9061620987350745659/32281802128991715328
        SC12=0.280703690307; SC7=0.2807037; SKN=43; SFL=24402
        SREPS=1513; SADM=29092719; SREJ=20485265; SMULT=5; SGRP=106; SLAB=14704790947
        SBLOB=718d7e8a3a3b937a519c3a04e14421032f34ec29ac99e9aa6d8bf9ab3b1e9ae2; SNSP=17; SNB=4 ;;
    K5lt) SD=sigma_catalogue/K5lt; SPRE=cosig3_s3_K5lt_f.cosig3; SKEYS=cosig3_s3_K5lt_f.keys.json; SPROB=s3_K5lt
        SPREP=s7_prep.py; SEVAL=s7_eval; SRN=K5lt_f
        SB=2747492282707/4123168604160; SSIG=1375676321453/4123168604160; SC12=0.333645420191; SC7=0.3336455
        SKN=45; SFL=28108; SREPS=1728; SADM=40527893; SREJ=16095211; SMULT=7; SGRP=108; SLAB=20294982806
        SZB=393705905561403120; SPC='"n_copies7": 630, "reps6": 1728, "labelled6": 868080'
        SBLOB=e73f987666b2a3582176116a94fbd6479ac9f3bc9a6a8cc4ccab2d32fa9b7ed1; SNSP=17; SNB=4 ;;
    K5_4minus) SD=sigma_catalogue/K5_4minus; SPRE=sigma4_c4_K5m_f.sigma4; SKEYS=sigma4_c4_K5m_f.keys.json
        SPROB=c4_K5m; SPREP=s7_prep.py; SEVAL=s7_eval; SRN=K5m4_f
        SB=24504905122391/30786325577728; SSIG=6281420455337/30786325577728; SC12=0.204032807991; SC7=0.2040329
        SKN=10; SFL=2243; SREPS=62; SADM=10058620; SREJ=54953092; SMULT=7; SGRP=56; SLAB=1741143596
        SZB=346363323387273360; SPC='"n_copies7": 105, "reps6": 62, "labelled6": 12068'
        SBLOB=7562c6add933e70414550ef2b83cf2149f169a8f9db412811d8fd6f29158bb02; SNSP=16; SNB=3 ;;
    K5eq) SD=sigma_catalogue/K5eq; SPRE=cosig3_s3_K5eq.cosig3; SKEYS=cosig3_s3_K5eq.keys.json; SPROB=s3_K5eq
        SPROG=s8; SPREP=s8_prepvm.py; SEVAL=s8_eval; SRN=K5eq_A
        SB=106600906627513/173173081374720; SSIG=66572174747207/173173081374720; SC12=0.384425652179; SC7=0.3844257
        SKN=45; SFL=29557; SREPS=1819; SADM=45962018; SREJ=13642974; SMULT=12; SGRP=108; SLAB=23083083434
        SZB=327925293541166160; SPC='"reps6": 1819, "labelled6": 918267'
        SBLOB=f598f7334ee035318cb29a094276225ea68d9ef9a169c41578aab401247a47ef; SNSP=17; SNB=4 ;;
    K5_3minus) SD=sigma_catalogue/K5_3minus; SPRE=cosig3_s3_K5m_f.cosig3; SKEYS=cosig3_s3_K5m_f.keys.json
        SPROB=s3_K5m; SPROG=s8; SPREP=s8_prepvm.py; SEVAL=s8_eval; SRN=K5m_f
        SB=1143076695549/1924145348608; SSIG=781068653059/1924145348608; SC12=0.405930172388; SC7=0.4059302
        SKN=46; SFL=32089; SREPS=1970; SADM=57357848; SREJ=7195112; SMULT=2; SGRP=118; SLAB=28854605654
        SZB=270167075410372560; SPC='"reps6": 1970, "labelled6": 990687'
        SBLOB=1aa762b5b9f82d0da2e9e46f8d20a45bcd77b6e45addbfa3c98c80272cb71b8e; SNSP=16; SNB=3 ;;
    C5) SD=sigma_catalogue/C5; SPRE=cosig3_s3_C5_f.cosig3; SKEYS=cosig3_s3_C5_f.keys.json; SPROB=s3_C5
        SPROG=s8; SPREP=s8_prepvm.py; SEVAL=s8_eval; SRN=C5_f
        SB=6422651703/8589934592; SSIG=2167282889/8589934592; SC12=0.252304935013; SC7=0.2523050
        SKN=40; SFL=13879; SREPS=835; SADM=5635379; SREJ=21725901; SMULT=1; SGRP=76; SLAB=2342442087
        SZB=226821172036026960; SPC='"reps6": 835, "labelled6": 371229'
        SBLOB=b3ae8f3197b93fbd22958ec32548e98c0d8c225b663f06270ac48eb5e1e0126d; SNSP=14; SNB=1 ;;
    K6_3) SD=sigma_catalogue/K6_3; SPRE=cosig3_s3_K6.cosig3; SKEYS=cosig3_s3_K6.keys.json; SPROB=s3_K6
        SPROG=s8; SPREP=s8_prepvm.py; SEVAL=s8_eval; SRN=K63_A
        SB=29749027005683/115448720916480; SSIG=85699693910797/115448720916480; SC12=0.742318262433; SC7=0.7423183
        SKN=48; SFL=35753; SREPS=2135; SADM=69958475; SREJ=1205; SMULT=11; SGRP=168; SLAB=34359509614
        SZB=310561986171245760; SPC='"reps6": 2135, "labelled6": 1048575'
        SBLOB=c6b8cfa21e13022311ab17c0923415454855c83356abc5a92e44e09cd5fec1d4; SNSP=17; SNB=4 ;;
    K6_4minus) SD=sigma_catalogue/K6_4minus; SPRE=sigma4_c4_K6m_f.sigma4; SKEYS=sigma4_c4_K6m_f.keys.json
        SPROB=c4_K6m; SPROG=s8; SPREP=s8_prepvm.py; SEVAL=s8_eval; SRN=K6m4_f
        SB=358317739722203/923589767331840; SSIG=565272027609637/923589767331840; SC12=0.612037993062; SC7=0.6120380
        SKN=12; SFL=6194; SREPS=154; SADM=160859017; SREJ=621687; SMULT=36; SGRP=149; SLAB=34244802014
        SZB=424828241254554480; SPC='"reps6": 154, "labelled6": 32752'
        SBLOB=ea9fb9fb75abe7a605fb6d41f68f3e2795d262304522961efee4b7f653f3ea78; SNSP=16; SNB=3 ;;
    K6_4) SD=sigma_catalogue/K6_4; SPRE=sigma4_c4_K6.sigma4; SKEYS=sigma4_c4_K6.keys.json; SPROB=c4_K6
        SPROG=s8; SPREP=s8_prepvm.py; SEVAL=s8_eval; SRN=K64_A
        SB=140873242475/549755813888; SSIG=408882571413/549755813888; SC12=0.743753064695; SC7=0.7437531
        SKN=12; SFL=6205; SREPS=155; SADM=162481445; SREJ=47835; SMULT=1; SGRP=164; SLAB=34352419335
        SZB=210085298235918000; SPC='"reps6": 155, "labelled6": 32767'
        SBLOB=349db8ec8680925e50ab7932008d993ce31d27beeac960eb294a834225c543aa; SNSP=11; SNB=1 ;;
    K7_4) SD=sigma_catalogue/K7_4; SPRE=sigma4_c4_K7.sigma4; SKEYS=sigma4_c4_K7.keys.json; SPROB=c4_K7
        SPROG=s8; SPREP=s8_prepvm.py; SEVAL=s8_eval; SRN=K74_A
        SB=9454851490593/76965813944320; SSIG=67510962453727/76965813944320; SC12=0.877155180904; SC7=0.8771552
        SKN=12; SFL=6206; SREPS=156; SADM=163577855; SREJ=1; SMULT=2; SGRP=174; SLAB=34359738367
        SZB=429451942557828960; SPC='"reps6": 156, "labelled6": 32768'
        SBLOB=76bc4012d7185cd1f6757098f31427a9e98adf369640db1d99e9322f433247f5; SNSP=14; SNB=4 ;;
  esac
}
SCNAMES="J4 K5lt K5_4minus K5eq K5_3minus C5 K6_3 K6_4minus K6_4 K7_4"
SCMASKS=0
sigcat_masks() {    # the spot-check graphs of verifier/sigcat_spot_values.txt as text and int64 files (once)
  [ "$SCMASKS" = 1 ] && return
  SCMASKS=1
  # shellcheck disable=SC2086
  runin sigcat sigcat_masks "$PY" -c '
import sys
import numpy as np
rows = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
for name in sys.argv[2:]:
    ms = [int(r[1]) for r in rows if r[0] == name]
    open("local/pc_%s.txt" % name, "w").write("".join("%d\n" % m for m in ms))
    np.array(ms, dtype=np.int64).tofile("local/pc_%s.bin" % name)
    print(name, len(ms), "masks written")' "$V/sigcat_spot_values.txt" $SCNAMES
  set --
  for n in $SCNAMES; do scinfo "$n"; set -- "$@" "$n $SNSP masks written"; done
  check sigcat_masks $? "$@"
}
SCSPOT=" "
sigcat_spot() {     # pure-Python exact values at the spot-check graphs, compared with the review's values (once)
  local n="$1"
  case "$SCSPOT" in *" $n "*) return ;; esac
  SCSPOT="$SCSPOT$n "
  scinfo "$n"
  runin sigcat "sigcat_pycheck_$n" "$PY" "$V/${SPROG}_pycheck.py" "$SPROB" "$C/$SD/$SPRE" "local/pc_$n.txt" \
    "local/pycheck_$n.json"
  check "sigcat_pycheck_$n" $?
  runin sigcat "sigcat_pyvals_$n" "$PY" -c '
import json, sys
from fractions import Fraction
name, spot, pyj, cert = sys.argv[1:5]
b = Fraction(json.load(open(cert))["bound"])
rows = [l.split() for l in open(spot) if l.strip() and not l.startswith("#")]
exp = {int(r[1]): Fraction(r[2]) for r in rows if r[0] == name}
got = {int(m): Fraction(v["slack"]) - b for m, v in json.load(open(pyj)).items()}
for m in exp:
    print("PY", m, "slack-b", got.get(m), "; review:", exp[m], "; equal:", got.get(m) == exp[m])
same = set(got) == set(exp) and all(got[m] == exp[m] for m in exp)
ge = all(v >= 0 for v in got.values())
print(name, len(got), "graphs; slack - b equal to the review value at every graph:", same, "; all >= b:", ge,
      "; at b:", sum(v == 0 for v in got.values()))
print("SPOT VALUES OK" if same and ge else "SPOT VALUES DIFFER")' "$n" "$V/sigcat_spot_values.txt" \
    "local/pycheck_$n.json" "$C/$SD/$SPRE.cert.json"
  check "sigcat_pyvals_$n" $? \
    "$n $SNSP graphs; slack - b equal to the review value at every graph: True ; all >= b: True ; at b: $SNB" \
    "SPOT VALUES OK"
}
SCBUILT=0
sigcat_build() {    # the exact evaluators (R7_SIG7: int64 and __int128 Gram entries; R7_SIG8: int64), once
  [ "$SCBUILT" = 1 ] && return
  SCBUILT=1
  mkdir -p "$W/sigcat/local"
  # shellcheck disable=SC2086
  run sigcat_build "$CC" $CFLAGS -o sigcat/s7_eval "$V/s7_eval.c"
  rc=$?
  # shellcheck disable=SC2086
  [ "$rc" = 0 ] && { ( cd "$W" && "$CC" $CFLAGS -o sigcat/s7_eval128 "$V/s7_eval128.c" ) >> "$W/logs/sigcat_build.log" 2>&1; rc=$?; }
  # shellcheck disable=SC2086
  [ "$rc" = 0 ] && { ( cd "$W" && "$CC" $CFLAGS -o sigcat/s8_eval "$V/s8_eval.c" ) >> "$W/logs/sigcat_build.log" 2>&1; rc=$?; }
  check sigcat_build $rc
}

if [ "$SIGCAT" = 1 ]; then
  # --- (m) Theorem sigmacat, fast part: key lists, reductions, exact values at the spot-check graphs, decimals -------
  set --; SCK=()
  for n in $SCNAMES; do
    scinfo "$n"; SCK+=("$C/$SD/$SKEYS")
    case "$SPRE" in *.cosig3) o=cosig3 ;; *) o=sigma4 ;; esac
    set -- "$@" "\"problem\": \"$SPROB\", \"objective\": \"$o\", \"keys\": $SKN, \"flags\": $SFL, \"bound\": \"$SB\", \"sigma_upper\": \"$SSIG\", \"PASS\": true"
  done
  run sigcat_keys "$PY" "$V/sk_check_keys.py" "${SCK[@]}"
  check sigcat_keys $? "ALL PASS: 10 key files checked" "$@"
  # the independent comparison of all ten key lists with the reviewer's own types and flags (R7_SIG8, addenda A1 for
  # the seven of that review and A2 for the three of R7_SIG7)
  set --; S3=()
  for n in $SCNAMES; do
    scinfo "$n"
    a=("$SRN {'hashes_equal': True, 'scalars_equal': True,"
       "'keys_equal': True, 'keys_vs_cert_json': True, 'n_keys': $SKN, 'total_flags': $SFL, 'counts_equal': True,")
    if [ "$SPROG" = s8 ]; then set -- "$@" "${a[@]}"; else S3+=("${a[@]}"); fi
  done
  runin sigcat sigcat_keys_indep7 "$PY" "$V/s8_keys_cmp.py" local/results_keys_cmp.json
  check sigcat_keys_indep7 $? "$@" \
    "controls (True = detected) {'flags_swapped': True, 'aut_changed': True, 'flag_replaced': True, 'rows_changed': True" \
    "ALL_PASS True" "!'PASS': False"
  runin sigcat sigcat_keys_indep3 "$PY" "$V/s8_keys_cmp_stage1.py" local/results_keys_cmp_stage1.json
  check sigcat_keys_indep3 $? "${S3[@]}" \
    "controls (True = detected) {'flags_swapped': True, 'aut_changed': True, 'flag_replaced': True, 'rows_changed': True" \
    "'K5lt_spec_replaced_by_K5eq': True, 'lam_2_to_1': True}" "ALL_PASS True" "!'PASS': False"
  mkdir -p "$W/sigcat/local"; rm -rf "$W/sigcat/local/neg_keys"
  runin sigcat sigcat_keys_neg "$PY" "$V/sk_neg_keys.py" local/neg_keys
  check sigcat_keys_neg $? "swap: exit 1 (FAIL as required)" "aut: exit 1 (FAIL as required)" \
    "forb: exit 1 (FAIL as required)" "rows: exit 1 (FAIL as required)" "NEGATIVE CONTROLS OK"
  runin sigcat sigcat_reduction "$PY" "$V/s7_reduction.py" results_reduction.json
  check sigcat_reduction $? '"r3_identities_graphs": 60' '"r3_num7_graphs": 40' '"r4_identities_graphs": 24' \
    '"r4_K5m_predicate_graphs": 900' '"planted_K4m_free_with_C5": 240' '"ALL_OK": true'
  runin sigcat sigcat8_predicates "$PY" "$V/s8_predicates.py" local/results_predicates8.json
  check sigcat8_predicates $? '"admissible": 284,' '"copies7": 315' '"admissible": 295,' '"copies7": 210' \
    '"admissible": 617,' '"admissible": 181,' '"copies7": 252' '"admissible": 169,' '"admissible": 323,' \
    '"admissible": 630,' '"copies7": 1' '"ALL_OK": true'
  sigcat_masks
  for n in $SCNAMES; do sigcat_spot "$n"; done
  # the exact bounds, 1 - b, the 12-digit ceilings of the reviews and the 7-digit ceilings of Table sigmacat
  # (upper bounds, rounded up)
  set --
  for n in $SCNAMES; do scinfo "$n"; set -- "$@" "$SD:$SPRE:$SKEYS:$SB:$SSIG:$SC12:$SC7"; done
  run sigcat_decimals "$PY" -c '
import json, os, sys
from decimal import Decimal, getcontext, ROUND_CEILING
from fractions import Fraction
getcontext().prec = 80
def ceil(x, q):
    return str((Decimal(x.numerator) / Decimal(x.denominator)).quantize(Decimal(q), rounding=ROUND_CEILING))
ok = True
for row in sys.argv[2:]:
    d, pre, keys, b, sig, c12, c7 = row.split(":")
    cj = json.load(open(os.path.join(sys.argv[1], d, pre + ".cert.json")))
    kj = json.load(open(os.path.join(sys.argv[1], d, keys)))
    B, S = Fraction(b), Fraction(sig)
    good = (Fraction(cj["bound"]) == B and Fraction(kj["bound"]) == B and Fraction(kj["sigma_upper"]) == S
            and 1 - B == S and ceil(1 - B, "1e-12") == c12 and ceil(1 - B, "1e-7") == c7)
    ok &= good
    print(d, "bound", B, "; sigma <= 1 - b =", 1 - B, "(paper:", sig + ") <", ceil(1 - B, "1e-12"), "(review:", c12 + "),",
          ceil(1 - B, "1e-7"), "(Table sigmacat:", c7 + ") ;", good)
print(len(sys.argv) - 2, "certificates;", "SIGCAT DECIMALS OK" if ok else "SIGCAT DECIMALS DIFFER")' "$C" "$@"
  check sigcat_decimals $? "10 certificates; SIGCAT DECIMALS OK"

  # --- (n) Remark j4limit: the obstruction witness (S1) and the exact dual point (S2) for sigma(J_4) ---------------
  rm -rf "$W/j4"; mkdir -p "$W/j4/local" "$W/j4/out"      # fresh caches (local/*.pkl) for every run
  runin j4 j4_s1 "$PY" "$V/jw_verify.py" "$C/sigma_catalogue/J4_limit/witness_Fanoc.json" --nullity --fullnull \
    --json out/verify_Fanoc.json
  check j4_s1 $? "24727 labelled 7-vertex masks, 93 classes" "E_phi f = 41/57, file target 41/57" \
    "|supp mu| = 87, outside phi: 0" "# keys 191, with kernel 187, sum kernel dims 1250, sum d(d+1)/2 6374, max |O| 156" \
    "[ok] V3_psd (0 failures)" "[ok] V3_entries_in_OxO" "[ok] V3_kernel" "[ok] V4_conditions (0 keys fail)" \
    "[ok] V5_delta_equals_file delta = 316891/14980607589" "6374 condition rows, rank mod p 86, nullity <= 7" \
    "[ok] N2_full_null_space_exact (all 7 basis vectors satisfy every condition exactly)" \
    "[ok] N3_vertex_optimal_exact_KKT" "[ok] N3_lp_optimum_equals_witness (LP delta = 316891/14980607589" \
    "VERDICT: PASS exit 0; delta = 316891/14980607589; fails = []" "![FAIL]"
  runin j4 j4_s1_controls "$PY" "$V/jw_neg.py"
  check j4_s1_controls $? "CONTROL OK  rv_moved " "CONTROL OK  rv_half " "CONTROL OK  rv_phi " \
    "CONTROL OK  producer_neg_phi " "CONTROL OK  mutate_edgeflip " "ALL CONTROLS AS EXPECTED" "!CONTROL BAD"
  runin j4 j4_s1_lift "$PY" "$V/jw_lift.py"
  check j4_s1_lift $? "two codings of the 9 blocks agree: True" \
    "lifting identities: 109374 entries checked on 33 graphs, 0 mismatches" "LIFT CHECK OK"
  runin j4 j4_s2 "$PY" "$V/jd_verify_dual.py" "$C/sigma_catalogue/J4_limit/dual_W2.json" --crosscheck \
    --relabel-check 12 --json out/verify_W2.json
  check j4_s2 $? "757 support graphs, den = 1298074214633706907132624082305024 (2^110: True)" \
    "[ok] D1_nums_positive" "[ok] D1_sum_equals_den" "[ok] D1_classes_distinct (757 distinct classes)" \
    "[ok] D1_all_J4_free (0 contain J_4)" "[ok] D2_V_equals_file" \
    "# V = 598379856892795294218114109955/831895706399775544732211478528 = 0.71929672468491; 41/57 - V = 1.5209291284e-06" \
    "# moments: 34 (block, canonical type) keys" "[ok] D3_crosscheck_multimodular_vs_bareiss (key (3, 5, 0), n = 130)" \
    "[ok] D3_all_moment_matrices_PSD (0 of 34 keys fail)" \
    "[ok] D4_relabel_lemma_spotcheck (1744 labelled (graph, block, type) keys, 0 mismatches)" "VERDICT: PASS; " "![FAIL]"
  # the window of Remark j4limit: weak duality between the J_4 certificate and the dual point, and the decimals
  run j4_window "$PY" -c '
import json, sys
from decimal import Decimal, getcontext, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction
getcontext().prec = 80
def dec(x, mode):
    return str((Decimal(x.numerator) / Decimal(x.denominator)).quantize(Decimal("1e-12"), rounding=mode))
V = Fraction(json.load(open(sys.argv[1]))["value"])
b = Fraction(json.load(open(sys.argv[2]))["bound"])
t = Fraction(41, 57)
print("V =", V, "; b of the J_4 certificate =", b)
print("b <= V (weak duality):", b <= V, "; V - b = %.4e" % float(V - b), "; 41/57 - V = %.4e" % float(t - V))
lo, hi = dec(1 - V, ROUND_FLOOR), dec(1 - b, ROUND_CEILING)
print("1 - V =", 1 - V, "> " + lo, "; 1 - b =", 1 - b, "< " + hi)
print("plain seven-vertex optimum for sigma(J_4) in (" + lo + ", " + hi + "]")
ok = b <= V < t and lo == "0.280703275315" and hi == "0.280703690307"
print("J4 WINDOW OK" if ok else "J4 WINDOW WRONG")' "$C/sigma_catalogue/J4_limit/dual_W2.json" \
    "$C/sigma_catalogue/J4/cosig3_s3_J4_r30.cosig3.cert.json"
  check j4_window $? "b <= V (weak duality): True ; V - b = 4.1499e-07 ; 41/57 - V = 1.5209e-06" \
    "plain seven-vertex optimum for sigma(J_4) in (0.280703275315, 0.280703690307]" "J4 WINDOW OK"
fi

# --- (p) Proposition lotup: upper bounds for lottery numbers (type-pattern certificates, Giraud host) -----------------
# The programs of the independent reviews R7_T3b (rv3b_*) and R7_P6 (r7p6_check.py), run in work/lotup/.
if [ "$LOTUP" = 1 ]; then
  mkdir -p "$W/lotup"
  runin lotup lotup_certs "$PY" "$V/rv3b_certs.py"
  check lotup_certs $? \
    "k7_r4p5_m5_den1000000.txt  [PASS]  sha256 5f0f88f5007a42535dfc17823cf3e1165749dbd9f27dcca8f1a0271fc0060d61 (sidecar match)" \
    "(k,r,p)=(7,4,5) m=5  types=38 {(4,): 5, (3, 1): 10, (2, 2): 10, (2, 1, 1): 12, (1, 1, 1, 1): 1}  blocks(y>0)=27  inner=None" \
    "Turan: multiset failures 0 ; brute force on 5x5 points: True (53130 5-sets)" "rows: 27/38 tight, min slack 0.0" \
    "c = 3467292853891919529774577/7750000000000000000000000 = 0.447392626309 ; stated c matches: True ; rounded-up 9 dp: 0.447392627 ; claimed decimal 0.447392627 ok: True" \
    "f8_k7_r4p6_x_quarter.txt  [PASS]  sha256 44abf7407bfa5dfc8b6d49e7aa8fadef33c8e563d39ca5c6325d752a54123403 (sidecar match)" \
    "(k,r,p)=(7,4,6) m=4  types=14 {(4,): 2, (2, 2): 6, (3, 1): 6}  blocks(y>0)=8  inner=[0, 1] q=5 c=3467292853891919529774577/7750000000000000000000000" \
    "Turan: multiset failures 0 ; brute force on 4x6 points: True (134596 6-sets)" "rows: 14/14 tight, min slack 0.0" \
    "c = 839651878561675758589323731/2976000000000000000000000000 = 0.282141088226 ; stated c matches: True ; rounded-up 9 dp: 0.282141089" \
    "ALL CERTIFICATES PASS: True" "negative controls k7_r4p5_m5_den1000000.txt" "negative controls f8_k7_r4p6_x_quarter.txt" \
    "!rejected: False" "![FAIL]"
  runin lotup lotup_derived "$PY" "$V/rv3b_derived.py"
  check lotup_derived $? "F8: 35*sum(y) = 107/384 (== 107/384: True) ; inner weight x0^4 + x1^4 = 1/128 (== 1/128: True)" \
    "c(7,4,6) <= 35 sum y + c745/128 with c745 = m5 cert: ... = 0.282141088226 ; == (107+3c)/384: True ; up 9dp 0.282141089" \
    "claimed 0.282141089 >= exact: True"
  runin lotup lotup_giraud "$PY" "$V/rv3b_giraud.py"
  check lotup_giraud $? "(T) d=2 n=8: Turan(5,4) by brute force: True" "(T) d=3 n=16: Turan(5,4) by brute force: True" \
    "(T) d=4 n=32: Turan(5,4) by brute force: True" \
    "(C) d=3: |E_cross| = 336 (formula 336) ; counts over ALL cross edges: [0] ; formula (N/2-2)(N/4-2) = 0" \
    "(C) d=4: |E_cross| = 6720 (formula 6720) ; counts over ALL cross edges: [12] ; formula (N/2-2)(N/4-2) = 12" \
    "(C) d=5: |E_cross| = 119040 (formula 119040) ; counts over ALL cross edges: [84] ; formula (N/2-2)(N/4-2) = 84" \
    "(I) d=3: internal count 6, C(N-4,2) = 6" "(I) d=4: internal count 66, C(N-4,2) = 66" \
    "(P) d=4: 8960 cross hyperedges (= |E_cross|*cnt/9 = 8960) ; max # through two cross edges = 2 ; alpha_cross = 1/6 ; internal (N-5)/C(N-4,2) = 1/6" \
    "(P) d=5: max # cross hyperedges through a (3,2) 5-set = 6 (hand bound N-2 = 30) ; alpha_cross = max(that,1)/cnt = 1/14" \
    "    leading ratio = 7/16" "    poly == direct formula at d=4,5: True" \
    "(E) d=4: cross 4-sets covered by gauge-trivial 6-sets: 6720 (= |E_cross| 6720: True) ; 5-sets with no block meeting them in >= 4 points: 0"
  runin lotup lotup_check "$PY" "$V/r7p6_check.py" "$ROOT"
  check lotup_check $? "PASS (b) condition (i): 126 types of size 5, failures 0" \
    "PASS (b) brute force on 5x5 blow-up: 53130 5-sets, failures 0" "PASS (b) C(7,4) sum y = c745 = file c" \
    "PASS (c) f8_k7_r4p6_x_quarter.txt condition (ii): 14 rows, all satisfied, 14 tight" \
    "PASS (c) f8_k7_r4p6_x_quarter.txt T and y = the paper's verbal description (arcs 1->2,2->3,3->4,4->1,1->3,2->4)" \
    "PASS (c) (107+3c745)/384 < 0.282141089 (value 0.282141088226369)" \
    "PASS (c) case analysis of the printed proof of (i) is exhaustive and correct" \
    "PASS (a) d=4 tiling: each of the 6720 cross edges has weight exactly 1 (good (3,3)-sets: 8960); internal C(N-4,2)=66" \
    "PASS table (7,4,6): C(k,4)/c = 124.051410 -> floor 124.05 (printed 124.05)" "SUMMARY: 69 checks, 0 FAIL" "!FAIL "
fi

# --- (k) to (n), long: the sigma raw scan, C10 over all labelled bases, the sigmacat raw scans, the J_4 controls ----
for x in $RAWX; do
  case "$x" in
    sigma)
      sigma_prepare
      NT=$THREADS; [ "$NT" -ge 1 ] 2>/dev/null || NT=1
      TOT=124928          # 122 representatives x 1024 high link patterns
      t0=$(date +%s); rcs=0; logs=""; tights=""; pids=""
      for i in $(seq 0 $((NT - 1))); do      # NT single-threaded processes on consecutive task ranges
        lo=$((TOT * i / NT)); hi=$((TOT * (i + 1) / NT))
        ( cd "$W/sigma" && ./rs_eval4s local/blob_sharp_sym.bin raw local/reps6_p5.bin "$lo" "$hi" "local/tight_raw_$i.u64" ) \
          > "$W/logs/sigma_raw_part$i.log" 2>&1 &
        pids="$pids $!"
        logs="$logs $W/logs/sigma_raw_part$i.log"; tights="$tights,local/tight_raw_$i.u64"
      done
      for pid in $pids; do wait "$pid" || rcs=1; done
      ELAPSED=$(( $(date +%s) - t0 ))
      # shellcheck disable=SC2086
      cat $logs > "$W/logs/sigma_raw.log"
      check sigma_raw $rcs "RANGE 0 " "RANGE $((TOT * (NT - 1) / NT)) $TOT"
      # shellcheck disable=SC2086
      runin sigma sigma_raw_analyze "$PY" "$V/rs_analyze.py" local/blob_sharp.bin.json 33/64 \
        "$C/sigma_K5_4/sharp_klp5.json" "${tights#,}" $logs
      check sigma_raw_analyze $? '"admissible": 86952880' '"rejected": 40973392' '"zero": 105' '"negative": 0' \
        '"min_slack": "33/64"' '"min_slack_equals_b": true' "\"next_slack_minus_b\": \"$MU0\"" '"groups": 113' \
        '"producer_groups_same_F": true' '"producer_groups_same_maxZ": true' '"tight_graphs": 105' \
        '"tight_distinct_labelled": 105' '"tight_all_giraud": true' '"tight_classes_hit": 10' ;;
    stab_c10)
      stab_prepare
      runin stab stab_c10_full "$PY" "$V/r_lg10_full.py"
      check stab_c10_full $? '"R_by_part_sizes": {"1,8": 9, "2,7": 2304, "3,6": 86016, "4,5": 516096, "0,9": 1}' \
        '"survivors_total": 15636107, "labelled_giraud_on_10_vertices_by_formula": 15636107, "total_equals_G10": true, "R_with_count_mismatch": [], "n_mismatch": 0, "sample_2000_explicit_and_struct_ok": true, "C10_FULL_PASS": true' ;;
    sigcat_K5eq|sigcat_K5_3minus|sigcat_C5|sigcat_K6_3|sigcat_K6_4minus|sigcat_K6_4|sigcat_K7_4)
      # Theorem sigmacat, the seven certificates of R7_SIG8: streaming decoding (evaluation table in blobs/, the
      # representatives in reps/, summary in prep/), spot values in C, exhaustive scan, comparison
      n="${x#sigcat_}"
      sigcat_build
      sigcat_masks
      sigcat_spot "$n"
      scinfo "$n"
      rm -f "$W/sigcat/local/reps6_$SPROB.npy" "$W/sigcat/local/reps6_${SPROB}_aut.npy"   # own representatives, anew
      runin sigcat "sigcat_prep_$n" "$PY" "$V/s8_prepvm.py" "$SPROB" "$C/$SD/$SPRE" "$n" .
      check "sigcat_prep_$n" $? "\"evaluator\": \"s8_eval\", \"keys\": $SKN, \"keys_match_cert\": true" \
        "\"z_bound_apriori_mine\": \"$SZB\", \"z_bound_apriori_cert\": \"$SZB\"" "$SPC" "\"blob_sha256\": \"$SBLOB\""
      runin sigcat "sigcat_each_$n" ./s8_eval "blobs/$n.blob" list "local/pc_$n.bin" "$THREADS" \
        "local/list_pc_$n.json" "local/each_$n.txt"
      check "sigcat_each_$n" $? "done list: admissible $SNSP rejected 0 bad 0"
      runin sigcat "sigcat_cmp_$n" "$PY" "$V/s8_cmp.py" "local/pycheck_$n.json" "local/each_$n.txt" \
        "$C/$SD/$SPRE.cert.json" "local/cmp_$n.json"
      check "sigcat_cmp_$n" $? "{\"n\": $SNSP, \"all_equal_and_ge_b\": true, \"n_at_min\": $SNB}"
      runin sigcat "sigcat_raw_$n" ./s8_eval "blobs/$n.blob" raw "reps/$SPROB.reps.bin" "$THREADS" "local/raw_$n.json"
      check "sigcat_raw_$n" $? "done raw: admissible $SADM rejected $SREJ bad 0"
      runin sigcat "sigcat_analyze_$n" "$PY" "$V/s8_analyze.py" "local/raw_$n.json" "$C/$SD/$SPRE.cert.json" \
        "prep/$n.prep.json" "local/results_raw_$n.json"
      check "sigcat_analyze_$n" $? '"bad": 0' "\"min_slack\": \"$SB\"" '"min_equals_claimed": true' \
        "\"min_multiplicity\": $SMULT" "\"n_groups_mine\": $SGRP" '"group_keys_identical": true' \
        '"group_maxZ_identical": true' '"raw_extensions_equal": true' "\"labelled_admissible_7vertex\": $SLAB" \
        "\"reps6\": $SREPS" "\"sigma_upper_exact\": \"$SSIG\"" "\"sigma_upper_ceil12\": \"$SC12\"" '"PASS": true' ;;
    sigcat_J4|sigcat_K5lt|sigcat_K5_4minus)      # Theorem sigmacat: decode, spot values in C, exhaustive scan
      n="${x#sigcat_}"
      sigcat_build
      sigcat_masks
      sigcat_spot "$n"
      scinfo "$n"
      rm -f "$W/sigcat/local/reps6_$SPROB.npy" "$W/sigcat/local/reps6_${SPROB}_aut.npy"   # own representatives, anew
      runin sigcat "sigcat_prep_$n" "$PY" "$V/$SPREP" "$SPROB" "$C/$SD/$SPRE" "local/$n"
      if [ "$SPREP" = s7_prep128.py ]; then
        check "sigcat_prep_$n" $? '"keys_match_cert": true, "flags_equal_json": true' '"max_diag_G_bits": 69' \
          '"z_bound_bits": 81, "S_bound_bits": 89' "\"reps6\": $SREPS" "\"blob_sha256\": \"$SBLOB\""
      else
        check "sigcat_prep_$n" $? '"keys_match_cert": true' "$SPC" \
          "\"z_bound_apriori_mine\": \"$SZB\", \"z_bound_apriori_cert\": \"$SZB\"" "\"blob_sha256\": \"$SBLOB\""
      fi
      runin sigcat "sigcat_each_$n" "./$SEVAL" "local/$n.blob" list "local/pc_$n.bin" "$THREADS" \
        "local/list_pc_$n.json" "local/each_$n.txt"
      check "sigcat_each_$n" $? "done list: admissible $SNSP rejected 0 bad 0"
      runin sigcat "sigcat_cmp_$n" "$PY" "$V/s7_cmp.py" "local/pycheck_$n.json" "local/each_$n.txt" \
        "$C/$SD/$SPRE.cert.json" "local/cmp_$n.json"
      check "sigcat_cmp_$n" $? "{\"n\": $SNSP, \"all_equal_and_ge_b\": true, \"n_at_min\": $SNB}"
      runin sigcat "sigcat_raw_$n" "./$SEVAL" "local/$n.blob" raw "local/$n.reps.bin" "$THREADS" "local/raw_$n.json"
      check "sigcat_raw_$n" $? "done raw: admissible $SADM rejected $SREJ bad 0"
      runin sigcat "sigcat_analyze_$n" "$PY" "$V/s7_analyze.py" "local/raw_$n.json" "$C/$SD/$SPRE.cert.json" \
        "local/$n.prep.json" "local/results_raw_$n.json"
      check "sigcat_analyze_$n" $? '"bad": 0' "\"min_slack\": \"$SB\"" '"min_equals_claimed": true' \
        "\"min_multiplicity\": $SMULT" "\"n_groups_mine\": $SGRP" '"group_keys_identical": true' \
        '"group_maxZ_identical": true' '"raw_extensions_equal": true' "\"labelled_admissible_7vertex\": $SLAB" \
        "\"reps6\": $SREPS" "\"sigma_upper_exact\": \"$SSIG\"" "\"sigma_upper_ceil12\": \"$SC12\"" '"PASS": true' ;;
    j4_controls)      # Remark j4limit: positive and negative controls of the dual-point checker
      mkdir -p "$W/j4/local"
      runin j4 j4_s2_controls "$PY" "$V/jd_controls.py"
      check j4_s2_controls $? "dead-flag graph search: " "CONTROL OK  producer ctrl_pos_Fanoc " \
        "CONTROL OK  producer ctrl_neg_lp " "CONTROL OK  rv phi law " "CONTROL OK  rv half mass on heaviest graph " \
        "ALL CONTROLS AS EXPECTED" "!CONTROL BAD" ;;
  esac
done

# --- (g) optional: first implementation's checker on the sharp N = 6 certificates -----------------------------------
if [ "$PRODN6" = 1 ]; then
  for x in "p5 1/4" "p6 1/10" "p5_l2 1/2"; do
    set -- $x
    run "producer_sharp6_$1" "$PY" "$ROOT/search/flagalg/verify_cert.py" "$C/n6_dual_points/cert_turan_r4_$1_N6_all_sharp.json"
    check "producer_sharp6_$1" $? "equal: True; verified >= claim: True"
  done
fi

# --- (o) optional: the producer's own checkers of the two files of Remark j4limit -----------------------------------
if [ "$PRODJ4" = 1 ]; then
  run producer_j4_s1 "$PY" "$ROOT/search/sigma_catalogue/verify_witness.py" \
    "$C/sigma_catalogue/J4_limit/witness_Fanoc.json"
  check producer_j4_s1 $? '"C4_delta_equals_file": true' "VERDICT: PASS.  delta = 316891/14980607589"
  run producer_j4_s2 "$PY" "$ROOT/search/sigma_catalogue/verify_dual3.py" "$C/sigma_catalogue/J4_limit/dual_W2.json"
  check producer_j4_s2 $? '{"PASS": true, "problem": "J4"' \
    '"support": 757, "value": "598379856892795294218114109955/831895706399775544732211478528"' \
    "PASS: no plain N = 7 flag-algebra certificate"
fi

say "# $NPASS passed, $NFAIL failed"
if [ "$NFAIL" = 0 ]; then say "ALL CHECKS PASSED"; exit 0; else say "SOME CHECKS FAILED"; exit 1; fi
