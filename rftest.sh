#!/bin/bash
# Automated tests for the RF cost layers (-CRF) and a few cross-cost regressions.
#
# Sections:
#   1. tree-vs-tree --ptreecost (RF, plus DC/DCE smoke)
#   2. naive --odtnaivecost on small networks
#   3. DP and BB on random tree-child networks (cross-validated against naive)
#   4. --bbtimestats output format
#   5. --printretusage output format
#   6. --general / non-tree-child guard
#   7. HC -CRF smoke
#   8. Display-tree sampling consistency
#
# Run from repo root after `make` (or `make test`). Override binary with
# SUPNET=./supnet256, etc. Set V=1 for verbose progress output.
#
# Exit code = number of failures.

set -u
SUPNET=${SUPNET:-./supnet}
V=${V:-0}

fails=0
checks=0

color_ok=$'\033[32m'
color_no=$'\033[31m'
color_off=$'\033[0m'
[ -t 1 ] || { color_ok=""; color_no=""; color_off=""; }

PASS() { checks=$((checks+1)); [ "$V" = "1" ] && echo "  ${color_ok}ok${color_off} $1"; }
FAIL() { fails=$((fails+1)); checks=$((checks+1)); echo "  ${color_no}FAIL${color_off} $*"; }

SECTION() { echo; echo "== $* =="; }

assert_eq() { # name expected actual
  if [ "$2" = "$3" ]; then PASS "$1"; else FAIL "$1: expected '$2', got '$3'"; fi
}

assert_match() { # name pattern actual
  if echo "$3" | grep -Eq "$2"; then PASS "$1"; else FAIL "$1: expected match '$2', got '$3'"; fi
}

assert_le() { # name a b  (a <= b for integers)
  if [ "$2" -le "$3" ]; then PASS "$1"; else FAIL "$1: $2 > $3"; fi
}

# -------- 1. tree-vs-tree --------
SECTION "1. ptreecost"

assert_eq "RF identical"     "0" "$($SUPNET --ptreecost -CRF -g '(a,(b,(c,d)))' -s '(a,(b,(c,d)))' 2>&1)"
assert_eq "RF one swap"      "2" "$($SUPNET --ptreecost -CRF -g '(a,(b,(c,d)))' -s '(a,((b,c),d))' 2>&1)"
assert_eq "RF balanced clash" "4" "$($SUPNET --ptreecost -CRF -g '((a,b),(c,d))' -s '((a,c),(b,d))' 2>&1)"
assert_eq "RF nontrivial"    "4" "$($SUPNET --ptreecost -CRF -g '(a,(b,(c,(d,e))))' -s '((a,b),((c,d),e))' 2>&1)"

# DC/DCE smoke (regression: don't break existing costs)
assert_eq "DC swap"          "1" "$($SUPNET --ptreecost -CDC  -g '(b,(a,c))' -s '(a,(b,c))' 2>&1)"
assert_eq "DCE swap"         "5" "$($SUPNET --ptreecost -CDCE -g '(b,(a,c))' -s '(a,(b,c))' 2>&1)"

# -------- 2. naive --------
SECTION "2. odtnaivecost"

assert_eq "naive tree-vs-tree" "2" "$($SUPNET -g '(a,(b,(c,d)))' -n '((a,b),(c,d))' -CRF --odtnaivecost 2>&1)"
assert_eq "naive 1-ret"        "2" "$($SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CRF --odtnaivecost 2>&1)"

# -------- 3. DP <= naive, BB == naive on random tree-child sweeps --------
SECTION "3. DP/BB vs naive (random tree-child sweeps)"

SEEDS="1 7 13 17 41 71 113 251 401 919 1009 4001"
GTS='(a,(b,(c,(d,(e,f)))))
((a,b),((c,d),(e,f)))
(((a,b),c),(d,(e,f)))
((a,b),(c,(d,(e,f))))
(((a,(b,c)),d),(e,f))
((((a,b),c),d),(e,f))'

for s in $SEEDS; do
  N=$($SUPNET -A6 -r1 -R3 --randseed $s --pnetworks)
  while IFS= read -r G; do
    [ -z "$G" ] && continue
    NAIVE=$($SUPNET -g "$G" -n "$N" -CRF --odtnaivecost 2>&1)
    DP=$($SUPNET    -g "$G" -n "$N" -CRF --DP          2>&1)
    BB=$($SUPNET    -g "$G" -n "$N" -CRF --BB          2>&1)
    case "$DP$BB$NAIVE" in *[!0-9]*)
      FAIL "non-integer in seed=$s G=$G N=$N: naive=$NAIVE DP=$DP BB=$BB"
      continue;;
    esac
    assert_eq "BB=naive seed=$s G=$G"    "$NAIVE" "$BB"
    assert_le "DP<=naive seed=$s G=$G"   "$DP"    "$NAIVE"
  done <<EOF
$GTS
EOF
done

# Higher reticulation count
for s in 7 31 91; do
  N=$($SUPNET -A8 -r1 -R5 --randseed $s --pnetworks)
  G='(a,(b,(c,(d,(e,(f,(g,h)))))))'
  NAIVE=$($SUPNET -g "$G" -n "$N" -CRF --odtnaivecost 2>&1)
  DP=$($SUPNET    -g "$G" -n "$N" -CRF --DP          2>&1)
  BB=$($SUPNET    -g "$G" -n "$N" -CRF --BB          2>&1)
  assert_eq "BB=naive 5ret seed=$s"  "$NAIVE" "$BB"
  assert_le "DP<=naive 5ret seed=$s" "$DP"    "$NAIVE"
done

# -------- 4. --bbtimestats output format --------
SECTION "4. --bbtimestats output format"

# Format: "<cost> <raw> <time> <minrt> <naivecnt> <naivetime> <dpcnt> <dptime>"
OUT=$($SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CRF --BB --bbtimestats 2>&1)
NCOL=$(echo "$OUT" | awk '{print NF}')
assert_eq "RF bbtimestats column count" "8" "$NCOL"
COST=$(echo "$OUT" | awk '{print $1}')
RAW=$(echo  "$OUT" | awk '{print $2}')
assert_eq "RF bbtimestats cost==raw" "$COST" "$RAW"

# DC: cost differs from raw (raw=DCE, cost=DC=DCE-2*lf-2)
OUT=$($SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CDC --BB --bbtimestats 2>&1)
DC=$(echo "$OUT" | awk '{print $1}')
DCE=$(echo "$OUT" | awk '{print $2}')
assert_eq "DC bbtimestats raw=DC+2*lf+2" "$((DC + 4*2 + 2))" "$DCE"

# -------- 5. --printretusage output format --------
SECTION "5. --printretusage output format"

OUT=$($SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CRF --DP --printretusage 2>&1)
NCOL=$(echo "$OUT" | awk '{print NF}')
assert_eq "RF DP printretusage cols" "2" "$NCOL"
assert_match "RF retusage format lft|rgh" '^[0-9]+ [0-9]+\|[0-9]+$' "$OUT"

OUT=$($SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CDC --DP --printretusage 2>&1)
assert_match "DC retusage format" '^[0-9-]+ [0-9]+\|[0-9]+$' "$OUT"

# Without flag: only the cost
OUT=$($SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CRF --DP 2>&1)
NCOL=$(echo "$OUT" | awk '{print NF}')
assert_eq "no flag => 1 column" "1" "$NCOL"

# -------- 6. tree-child guard --------
SECTION "6. non-tree-child guard"

# Network (((a)#X,(b)#Y),((#X,#Y),(c,d))): two nodes with all-ret children.
NTC='(((a)#X,(b)#Y),((#X,#Y),(c,d)))'
RAW=$($SUPNET -n "$NTC" --detectclass 2>&1)
TC_FLAG=$(echo "$RAW" | awk '{print $2}')
assert_eq "test fixture is non-tree-child" "0" "$TC_FLAG"

# DP must bail out (and exit nonzero) on non-tree-child input
$SUPNET -g '(a,(b,(c,d)))' -n "$NTC" -CRF --DP --general >/dev/null 2>/tmp/rftest_err.$$
RC=$?
if [ $RC -ne 0 ]; then PASS "non-tree-child --DP exits nonzero"; else FAIL "non-tree-child --DP returned 0"; fi
assert_match "non-tree-child error message" 'tree-child' "$(cat /tmp/rftest_err.$$ 2>/dev/null)"
rm -f /tmp/rftest_err.$$

# DC DP regression: still works on a regular tree-child network.
$SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CDC --DP >/dev/null 2>&1
RC=$?
if [ $RC -eq 0 ]; then PASS "DC DP regression on tree-child net"; else FAIL "DC DP failed (rc=$RC)"; fi

# -------- 7. HC -CRF smoke --------
SECTION "7. HC -CRF smoke"

OUT=$($SUPNET -A4 -r1 --randseed 1 -CRF --HC --noodtfiles 2>&1 | tail -1)
assert_match "HC RF summary" 'Cost:[0-9]+.*Class:TreeChild' "$OUT"

OUT=$($SUPNET -g '(a,(b,(c,d)))' -n '((#1,d),(((c)#1,b),a))' -CRF --HC --noodtfiles 2>&1 | tail -1)
assert_match "HC RF on fixed input" 'Cost:[0-9]+' "$OUT"

# -------- 8. Display-tree sampling consistency --------
SECTION "8. display-tree sampling"

NSAMP='((a)#X,((b,(c,(d,e))),(#X,f)))'
GSAMP='(a,(b,(c,(d,(e,f)))));((a,b),((c,d),(e,f)));(((a,b),c),(d,(e,f)))'
NAIVE=$($SUPNET -g "$GSAMP" -n "$NSAMP" -CRF --odtnaivecost 2>&1)
for sr in 0.3 0.5 1.0; do
  for seed in 1 7 13; do
    OUT=$($SUPNET -g "$GSAMP" -n "$NSAMP" -CRF --odtnaivecost \
            --displaytreesampling $sr --randseed $seed 2>&1)
    assert_le "sampling rate=$sr seed=$seed >= naive" "$NAIVE" "$OUT"
  done
done

# -------- summary --------
echo
if [ $fails -eq 0 ]; then
  echo "${color_ok}rftest.sh: $checks checks, 0 fails${color_off}"
else
  echo "${color_no}rftest.sh: $checks checks, $fails fails${color_off}"
fi
exit $fails
