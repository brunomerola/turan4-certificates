#!/usr/bin/env bash
# Run the independent checks of this package from its root directory (see README.md and verifier/EXPECTED.txt).
#
#   verifier/verify_all.sh [--fast] [--raw NAMES] [--producer-n6] [--threads N] [--work DIR]
#
#   --fast          every check except the raw scans (the default when no option is given); about 30 min on 1 core
#   --raw NAMES     the full exhaustive raw scan (a) for the named 4-graph certificates, comma-separated, from
#                   K5_4, K6_4, K7_4, K5_4minus, K6_4minus, or "all". K5_4minus takes minutes; each of the
#                   others takes 1 to 3.5 hours on 2 threads (about 5 minutes on 28 threads)
#   --producer-n6   also run the FIRST implementation's checker search/flagalg/verify_cert.py on the three sharp
#                   N = 6 certificates (this is not an independent check; see README.md)
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
FAST=0; RAW=""; PRODN6=0; THREADS=2; WORK="work"
while [ $# -gt 0 ]; do
  case "$1" in
    --fast) FAST=1 ;;
    --raw) RAW="${2:?--raw needs a list}"; shift ;;
    --producer-n6) PRODN6=1 ;;
    --threads) THREADS="${2:?--threads needs a number}"; shift ;;
    --work) WORK="${2:?--work needs a directory}"; shift ;;
    -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "unknown option: $1 (see --help)"; exit 2 ;;
  esac
  shift
done
if [ "$FAST" = 0 ] && [ -z "$RAW" ] && [ "$PRODN6" = 0 ]; then FAST=1; fi
[ "$RAW" = all ] && RAW="K5_4,K6_4,K7_4,K5_4minus,K6_4minus"
RAW="${RAW//,/ }"
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

# Parameters of the five 4-graph certificates (paper: Theorem main, Table counts):
# directory, file prefix, p, lambda, 6-vertex representatives, classes on 7 vertices, labelled admissible graphs on
# 7 vertices, admissible raw extensions, exact bound, certificate argmin (colex), raw-scan NMIN numerator and
# denominator 5040 M^2, raw-scan checksum (sumZ, sume).
cinfo() {
  case "$1" in
    K5_4) D=K5_4; PRE=lpcg_p5_conv_full; P=5; LAM=1; NREP=122; NCLS=3908438; NLAB=19199206747; NRAW=86952880
          FRAC=35604499940047/115448720916480; ARG=27303369217; NMIN=1709015997122256; DEN=5541538603991040
          SUMZ=3261680919295407098048; SUME=1670012489 ;;
    K6_4) D=K6_4; PRE=lpcg_p6_full; P=6; LAM=1; NREP=155; NCLS=7011184; NLAB=34352419335; NRAW=162481445
          FRAC=614350111275/4398046511104; ARG=34359738367; NMIN=3096324560826000; DEN=22166154415964160
          SUMZ=145576489784956352362440; SUME=2851672393 ;;
    K7_4) D=K7_4; PRE=lpcg_p7_full; P=7; LAM=1; NREP=156; NCLS=7013319; NLAB=34359738367; NRAW=163577855
          FRAC=7476698908057/115448720916480; ARG=33940849153; NMIN=358881547586736; DEN=5541538603991040
          SUMZ=97827306173207768161920; SUME=2862612480 ;;
    K5_4minus) D=K5_4minus; PRE=lpcg_h44e; P=5; LAM=2; NREP=62; NCLS=360127; NLAB=1741143596; NRAW=10058620
          FRAC=2430536617277/4398046511104; ARG=34359738367; NMIN=12249904551076080; DEN=22166154415964160
          SUMZ=5935335888043916906616; SUME=230039592 ;;
    K6_4minus) D=K6_4minus; PRE=lpcg_k6m_f; P=6; LAM=2; NREP=154; NCLS=6986573; NLAB=34244802014; NRAW=160859017
          FRAC=964431985683/4398046511104; ARG=31837153552; NMIN=4860737207842320; DEN=22166154415964160
          SUMZ=537778008513653740860644; SUME=2833905684 ;;
    *) echo "unknown certificate name: $1"; exit 2 ;;
  esac
  if [ "$LAM" = 1 ]; then SUF=""; else SUF="_l$LAM"; fi
}

say "# verify_all.sh  $(date -u +%FT%TZ)  root=$ROOT  work=$W  threads=$THREADS  python=$("$PY" -c 'import sys; print(sys.version.split()[0])' 2>&1)"
say "# options: fast=$FAST raw='${RAW}' producer_n6=$PRODN6"

# --- 0. integrity and build -----------------------------------------------------------------------------------------
if command -v sha256sum > /dev/null; then SHA="sha256sum"; else SHA="shasum -a 256"; fi
t0=$(date +%s); (cd "$ROOT" && $SHA -c SHA256SUMS) > "$W/logs/integrity.log" 2>&1; rc=$?; ELAPSED=$(( $(date +%s) - t0 ))
check integrity $rc "!FAILED"
run npz_hashes "$PY" -c '
import hashlib, json, sys, glob, os
ok = True
for j in sorted(glob.glob(os.path.join(sys.argv[1], "*", "*.cert.json"))):
    d = json.load(open(j))
    h = hashlib.sha256(open(j[:-5] + ".npz", "rb").read()).hexdigest()
    print(os.path.basename(j), "cert_npz_sha256", d["cert_npz_sha256"] == h)
    ok &= d["cert_npz_sha256"] == h
print("ALL NPZ HASHES MATCH" if ok else "NPZ HASH MISMATCH")' "$C"
check npz_hashes $? "ALL NPZ HASHES MATCH"
NEED_EVAL=0
{ [ "$FAST" = 1 ] || [ -n "$RAW" ]; } && NEED_EVAL=1
if [ "$NEED_EVAL" = 1 ]; then
  # shellcheck disable=SC2086
  run build "$CC" $CFLAGS -o rv_eval "$V/rv_eval.c"
  check build $?
fi

# --- (a) the five 4-graph certificates ------------------------------------------------------------------------------
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

ALL4="K5_4 K6_4 K7_4 K5_4minus K6_4minus"
if [ "$FAST" = 1 ]; then
  for n in $ALL4; do
    prepare "$n"; cinfo "$n"
    run "pycheck_cert_$n" env RV_LAM="$LAM" "$PY" "$V/rv_pycheck.py" "$C/$D/$PRE" "$ARG"
    check "pycheck_cert_$n" $? "$ARG $FRAC " "== cert bound: True"
  done
fi
for n in $RAW; do
  prepare "$n"; cinfo "$n"
  run "raw_$n" env RV_LAM="$LAM" OMP_NUM_THREADS="$THREADS" ./rv_eval "blob_$n.bin" raw "reps_p${P}${SUF}.bin"
  check "raw_$n" $? "TOTAL admissible=$NRAW " "NMIN $NMIN DEN $DEN " "SUMS sumZ $SUMZ sume $SUME"
  run "compare_raw_$n" env RV_LAM="$LAM" "$PY" "$V/rv_compare.py" "$C/$D/$PRE.cert.json" "logs/raw_$n.log" --raw
  check "compare_raw_$n" $? "rv bound $FRAC = " "EQUAL: True" "same set: True; maxZ equal on all: True"
  RARG=$(awk '$1 == "NMIN" {print $8}' "$W/logs/raw_$n.log")
  run "pycheck_raw_$n" env RV_LAM="$LAM" "$PY" "$V/rv_pycheck.py" "$C/$D/$PRE" "${RARG:-0}"
  check "pycheck_raw_$n" $? "${RARG:-0} $FRAC " "== cert bound: True"
done

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
fi

# --- (g) optional: first implementation's checker on the sharp N = 6 certificates -----------------------------------
if [ "$PRODN6" = 1 ]; then
  for x in "p5 1/4" "p6 1/10" "p5_l2 1/2"; do
    set -- $x
    run "producer_sharp6_$1" "$PY" "$ROOT/search/flagalg/verify_cert.py" "$C/n6_dual_points/cert_turan_r4_$1_N6_all_sharp.json"
    check "producer_sharp6_$1" $? "equal: True; verified >= claim: True"
  done
fi

say "# $NPASS passed, $NFAIL failed"
if [ "$NFAIL" = 0 ]; then say "ALL CHECKS PASSED"; exit 0; else say "SOME CHECKS FAILED"; exit 1; fi
