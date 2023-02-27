#! /bin/bash
cluster=$1

convert results/$cluster/results/figures/labeled.png -crop 500x500+110+15 results/$cluster/results/figures/labeled_no_border.png
