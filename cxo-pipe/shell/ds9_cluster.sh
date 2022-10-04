#! /bin/bash
shopt -s expand_aliases
mapfile=$1

ds9 $mapfile -scale log -cmap b &
