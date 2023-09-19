#! /bin/usr/env bash

data='../../../data/210108_A00405_0331_BHNW52DSXY'
res='..'
sub='20210517_sheridan_geo'
sums='geo_md5sums.txt'
meta='geo_metadata.xlsx'

set -o errexit -o pipefail -o nounset -x

mkdir -p "$sub"

# geo metadata
ln -sr "$meta" "$sub"

# cellranger matrices
tmp=$(mktemp tmp.XXXXX)

arr=("$res/M?/" "$res/A?/")

for dir in ${arr[@]}
do
    nm="$(basename $dir)"
    mat_dir='outs/filtered_feature_bc_matrix'

    for file in 'barcodes.tsv.gz' 'features.tsv.gz' 'matrix.mtx.gz'
    do
        mat="$sub/${nm}_$file"

        ln -sr "$dir/$mat_dir/$file" "$mat"

        md5sum "$mat" \
            >> "$tmp"
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
for dir in ${arr[@]};
do
    fq="$(basename $dir)"
    fq=$data/$fq*.fastq.gz

    ln -sr $fq "$sub"
done

# format md5sums
cat $data/md5sums.txt \
    >> $tmp

cat $tmp \
    | awk -v OFS="  " '{gsub("^.*/", "", $2); print}' \
    | sort -k2,2 \
    > $sums

rm $tmp


