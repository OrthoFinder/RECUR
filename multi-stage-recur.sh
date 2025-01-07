#!/bin/bash

# if [ -z "$1" ] || [ -z "$2" ]; then
#     echo "Usage: $0 <NALIGN> <SEED0>"
#     echo "Please provide the number of alignments (NALIGN) and the initial seed (SEED0) as arguments."
#     exit 1
# fi

# NALIGN=$1
# SEED0=$2

NALIGN=10 # Number of alignments for each run
SEED0=10 # Starting seed number

MAXNALIGN=1000
if [ $((MAXNALIGN % NALIGN)) -ne 0 ]; then
    echo "Error: NALIGN must divide MAXNALIGN ($MAXNALIGN) evenly."
    exit 1
fi


SEED1=$((SEED0 + (MAXNALIGN / NALIGN) - 1)) 

# -ds stands for disk saving
# -ms stands for multi stage
# -ds can be used alone, however, -ms cannot, it needs to be used together with -ds
recur -f ExampleData -st AA --outgroups ExampleData --num-alignments $NALIGN -ds -ms --seed $SEED0

for i in $(seq $((SEED0 + 1)) $SEED1); do
    echo
    custom_command="recur -f ExampleData -st AA --outgroups ExampleData --num-alignments $NALIGN -ds -ms -rs 3 --seed $i"
    echo ">>>>> Run command: $custom_command"

    custom_result=$($custom_command 2>&1)
    if [ $? -ne 0 ]; then
        echo "Error executing command: $custom_command"
        echo "Error output:"
        echo "$custom_result"
        continue 
    fi

    echo
    echo "Command output:"
    echo "$custom_result"
    echo
done

recur -f ExampleData -cr

