#!/bin/bash

qsub -cwd -b y -q bhklab -e sge_out -o sge_out -N "cluster" -t 1-2016 -tc 15 "module load R; Rscript batch.cluster.all.R"
