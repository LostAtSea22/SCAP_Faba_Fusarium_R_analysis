#!/bin/bash
set -e  # stop if any command fails

Rscript 01_preprocess.R
Rscript 02_analysis.R
Rscript 03_plots.R
