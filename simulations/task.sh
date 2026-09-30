#!/bin/bash
NSEEDS=100
idx=$((HQ_TASK_ID - 1))
param_idx=$(( idx / NSEEDS + 1 ))
seed=$(( idx % NSEEDS + 1 ))
PARAMS=$(sed -n "${param_idx}p" params.txt)

Rscript --no-save --slave sim_areal.R $seed $PARAMS
