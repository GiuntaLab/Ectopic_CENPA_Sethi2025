#!/bin/sh

#AS-HOR-vs-RPE1v1.1.bed is the centromere annotation file provided by HumAS-HMMER_for_AnVIL
bedtools intersect -a $1 -b AS-HOR-vs-RPE1v1.1.bed -wa -wb > intersect_in_out_centromere/$1\_AS_HOR.bed
bedtools intersect -v -a $1 -b AS-HOR-vs-RPE1v1.1.bed -wa -wb > intersect_in_out_centromere/$1\_out_AS_HOR.bed
