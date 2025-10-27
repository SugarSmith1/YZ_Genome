#!/bin/bash

# Pore-C heatmap for all chromosomes
haphic plot YZ08169.chr.sorted.agp ./contact_matrix.pkl --threads 8 --min_len 10 --separate_plots

# Pore-C heatmap for chromosome 1
haphic plot chrs.agp contact_matrix.pkl