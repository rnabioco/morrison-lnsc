#! /bin/usr/env bash

data=('210416_A00405_0383_AHY2N3DSXY' '210716_A00405_0430_BHHFF2DSX2')
dat_dir='../../data'
res_dir=('../2021-04-16' '../2021-07-16')
sub='20230919_sheridan_geo'
sums='geo_md5sums.txt'
meta='geo_metadata.xlsx'

set -o errexit -o pipefail -o nounset -x

mkdir -p "$sub"

# geo metadata
ln -sr "$meta" "$sub"

# cellranger matrices
tmp=$(mktemp tmp.XXXXX)

nms=('AF1' 'AF2' 'M1' 'M2')

for res in ${res_dir[@]}
do
    for nm in ${nms[@]}
    do
        run=$(basename $res)
        dir="$res/$nm"
        mat_dir="$res/$nm/outs/filtered_feature_bc_matrix"
    
        for file in 'barcodes.tsv.gz' 'features.tsv.gz' 'matrix.mtx.gz'
        do
            mat="$sub/${run}_${nm}_$file"
    
            ln -sr "$mat_dir/$file" "$mat"
    
            md5sum "$mat" \
                >> "$tmp"
        done
    done
done

# Seurat metadata
for file in *_count_matrix.tsv.gz *_metadata.tsv.gz
do
    ln -sr "$file" "$sub"

    md5sum "$file" \
        >> "$tmp"
done

# fastqs
for dat in ${data[@]}
do
    fq=$dat_dir/$dat/*.fastq.gz

    ln -sr $fq "$sub"

    cat "$dat_dir/$dat/md5sums.txt" \
        | sort -k2,2 \
        | awk '$2 !~ ".csv$"' \
        >> $tmp
done

# format md5sums
cat $tmp \
    | awk -v OFS="  " '{gsub("^.*/", "", $2); print}' \
    > $sums

rm $tmp

