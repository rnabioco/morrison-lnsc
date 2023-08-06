#! /bin/usr/env bash

set -o errexit -o nounset -o pipefail -x

target_dir="$HOME/Projects/morrison-scRNA-seq/results"
dates=('2021-04-16' '2021-07-16')

for d in ${dates[@]}
do
    path="$target_dir/$d"
    paths=($(find "$path" -regextype 'posix-awk' -regex "$path/(A|AF|M)[0-9]+"))

    for p in ${paths[@]}
    do
        sam=$(basename "$p")
        mat='outs/filtered_feature_bc_matrix'

        mkdir -p "$d/$sam/$mat"

        rsync -vr "$p/$mat/" "$d/$sam/$mat"
    done
done

