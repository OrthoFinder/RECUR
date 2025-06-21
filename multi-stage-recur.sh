#!/usr/bin/env bash
set -euo pipefail

# HELP
show_help() {
cat <<'EOF'
Usage: ./multi-stage-recur.sh [OPTIONS]

  -d DATASET      dataset directory / prefix        (default: ExampleData)
  -n TOTAL_ALIGN  total Monte-Carlo simulations     (overrides Log.txt)
  -b ALIGN_BATCH  alignments per batch              (default: 1000)
  -s SEED0        starting seed                     (default: 8)
  -m PV_METHOD    p-value adjust method             (bonferroni, fdr_bh,…; overrides Log.txt)
  -e              enable MC error control           (adds -mce)
  -p              print p-value statistics          (adds -ps)
  -h              show this help and exit
EOF
exit 1
}

# DEFAULTS
TOTAL_ALIGN=""
TOTAL_REC=""
ALIGN_BATCH=1000
SEED0=8
DATASET="ExampleData"
PV_METHOD=""
USE_MCE=false
USE_PV_STATS=false

# OPTION PARSING
while getopts ":n:b:s:d:m:ep h" opt; do
  case "$opt" in
    n) TOTAL_ALIGN=$OPTARG ;;
    b) ALIGN_BATCH=$OPTARG ;;
    s) SEED0=$OPTARG ;;
    d) DATASET=$OPTARG ;;
    m) PV_METHOD=$OPTARG ;;
    e) USE_MCE=true ;;
    p) USE_PV_STATS=true ;;
    h) show_help ;;
    *) echo "Unknown flag: -$OPTARG" >&2; show_help ;;
  esac
done

echo ""
# A) Real_Phylogeny step (runs once, if needed)
real_dir=$(ls -d "$DATASET"/*.recur/Real_Phylogeny 2>/dev/null | head -1 || true)
if [[ -z $real_dir || ! -e $real_dir/*.treefile ]]; then
  echo "==== Running Real Phylogeny step ===="
  cmd=(recur -f "$DATASET" -st AA --outgroups "$DATASET" -ds -ms)
  [[ -n $PV_METHOD ]] && cmd+=(-pam "$PV_METHOD")
  $USE_MCE && cmd+=(-mce)
  echo ">>>>> ${cmd[*]}"
  "${cmd[@]}"
else
  echo "✓ Real Phylogeny already present — skipping"
fi
echo

# B) Read counts / PV_METHOD from Log.txt if missing
if [[ -z $TOTAL_ALIGN || -z $TOTAL_REC || -z $PV_METHOD ]]; then
  LOGFILE=$(ls -t "$DATASET"/*.recur/Log.txt 2>/dev/null | head -1 || true)
  [[ -z $LOGFILE ]] && { echo "No Log.txt found and counts not supplied."; exit 1; }

  read -r LOG_REC LOG_ALIGN LOG_PVM < <(
    grep -m1 -E 'MCS Parameters[[:space:]]*\|' "$LOGFILE" |
    awk 'BEGIN{FS="[[:space:]]*\\|[[:space:]]*"}
         {for(i=1;i<=NF;i++){split($i,a,"=");
           if(a[1]=="num_rec_tests")nr=a[2];
           if(a[1]=="num_mc_sims")  mc=a[2];
           if(a[1]=="pval_adj_method")pvm=a[2];}
          if(mc)print nr,mc,pvm}')
  [[ -z $TOTAL_REC   ]] && TOTAL_REC=$LOG_REC
  [[ -z $TOTAL_ALIGN ]] && TOTAL_ALIGN=$LOG_ALIGN
  [[ -z $PV_METHOD && -n $LOG_PVM ]] && PV_METHOD=$LOG_PVM
fi

# C) Sanity checks
for v in TOTAL_ALIGN TOTAL_REC ALIGN_BATCH SEED0; do
  [[ ${!v} =~ ^[0-9]+$ ]] || { echo "$v must be integer." >&2; exit 1; }
done
(( TOTAL_ALIGN > 0 )) || { echo "TOTAL_ALIGN must be >0." >&2; exit 1; }

# D) Derived values & summary banner
if (( TOTAL_ALIGN <= ALIGN_BATCH )); then
  NBATCH=0
  RES_ALIGN=$TOTAL_ALIGN
else
  NBATCH=$(( TOTAL_ALIGN / ALIGN_BATCH ))
  RES_ALIGN=$(( TOTAL_ALIGN % ALIGN_BATCH ))
fi

printf "\n====================  RUN SUMMARY  ====================\n"
printf "Dataset (-d)                    : %s\n"   "$DATASET"
printf "Recurrence tests                : %'d\n"  "$TOTAL_REC"
printf "Monte-Carlo simulations (-n)    : %'d\n"  "$TOTAL_ALIGN"
printf "Batch size (-b)                 : %'d\n"  "$ALIGN_BATCH"
printf "Number of full batches          : %'d\n"  "$NBATCH"
printf "Remainder simulations           : %'d\n"  "$RES_ALIGN"
printf "Starting seed (-s)              : %'d\n"  "$SEED0"
printf "P-value method (-m)             : %s\n"    "${PV_METHOD:-<none>}"
printf "MC error control (-e)           : %s\n"    "$($USE_MCE && echo enabled || echo disabled)"
printf "Print PV stats (-p)             : %s\n"    "$($USE_PV_STATS && echo enabled || echo disabled)"
printf "========================================================\n\n"

# E) One-shot OR multi-stage execution
if (( NBATCH == 0 )); then
  echo "==== One-shot run (TOTAL_ALIGN ≤ ALIGN_BATCH) ===="
  cmd=(recur -f "$DATASET" -st AA --outgroups "$DATASET" -ds -ms
       --num-alignments "$TOTAL_ALIGN" --mcs-seed "$SEED0")
  [[ -n $PV_METHOD ]] && cmd+=(-pam "$PV_METHOD")
  $USE_MCE && cmd+=(-mce)
  echo ">>>>> ${cmd[*]}"
  "${cmd[@]}"

else
  # base command for stage-3 batches
  base_cmd=(recur -f "$DATASET" -st AA --outgroups "$DATASET" -rs 3 -ds -ms)
  [[ -n $PV_METHOD ]] && base_cmd+=(-pam "$PV_METHOD")
  $USE_MCE && base_cmd+=(-mce)

  # full batches
  for ((i=0;i<NBATCH;i++)); do
    SEED=$(( SEED0 + i ))
    cmd=("${base_cmd[@]}" --num-alignments "$ALIGN_BATCH" --mcs-seed "$SEED")
    echo ">>>>> Batch $i: ${cmd[*]}"
    "${cmd[@]}"
  done

  # remainder
  if (( RES_ALIGN )); then
    SEED=$(( SEED0 + NBATCH ))
    cmd=("${base_cmd[@]}" --num-alignments "$RES_ALIGN" --mcs-seed "$SEED")
    echo ">>>>> Remainder: ${cmd[*]}"
    "${cmd[@]}"
  fi
fi

# F) Final combined recurrence list
echo -e "\n==== Generating combined recurrence list (-cr) ===="
final=(recur -f "$DATASET" -cr)
[[ -n $PV_METHOD ]] && final+=(-pam "$PV_METHOD")
$USE_PV_STATS && final+=(-ps)
echo ">>>>> ${final[*]}"
"${final[@]}"
