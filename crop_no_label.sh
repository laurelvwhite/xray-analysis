#! /bin/bash
cluster=$1

convert ../results/$cluster/results/figures/unlabeled.png -crop 500x500+110+15 ../results/$cluster/results/figures/unlabeled_no_border.png
