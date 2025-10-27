#!/bin/bash

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Creation Time: 2024-07-25 12:18


# conda activate jcvi

# mv Erianthus_rufipilus.cds Eru.cds
# python3 -m jcvi.formats.gff bed --type=mRNA --key=ID Erianthus_rufipilus.gff3 -o Eru.bed
python3 -m jcvi.formats.gff bed --type=mRNA --key=ID YZ081609.gff3 -o YZ081609.all.bed
grep "^chr" YZ081609.all.bed > YZ081609.bed

# python -m jcvi.compara.catalog ortholog YZ081609 Eru --no_strip_names