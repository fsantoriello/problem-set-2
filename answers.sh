#! /usr/bin/env bash

datasets="$HOME/Desktop/School_Stuff/Courses/GENOMICS/data-sets"

#Use BEDtools intersect to identify the size of the largest overlap 
#between CTCF and H3K4me3 locations.

CTCF="$datasets/bed/encode.tfbs.chr22.bed.gz"

H3K4me3="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"

answer_1=$(intersectBed -a $CTCF -b $H3K4me3 \
    | awk 'BEGIN {OFS="\t"} ($4=="CTCF") {print $0, $3-$2}' \
    | sort -k5nr \
    | cut -f 5 \
    | head -n 1)

echo "answer-1: $answer_1"

#Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
#19,000,500 on chr22 of 'hg19' genome build. Report the GC content as a
#fraction (e.g., 0.50).

chr22="$datasets/fasta/hg19.chr22.fa"

echo -e "chr22\t19000000\t19000500\n" > tmp.bed

answer_2=$(nucBed -fi $chr22 -bed tmp.bed \
    | grep -v '^#' \
    | cut -f 5)

rm tmp.bed

echo "answer-2: $answer_2"

#Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
#interval) that has the largest mean signal in 'ctcf.hela.chr22.bg.gz'.

CTCF="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

TFBS="$datasets/bed/encode.tfbs.chr22.bed.gz"

answer_3=$(mapBed -a $TFBS -b $CTCF -c 4 -o mean \
    | awk '($4=="CTCF")' \
    | sort -k5nr \
    | awk 'BEGIN {OFS="\t"} (NR==1) {print $0, $3-$2}' \
    | cut -f 6)

echo "answer-3: $answer_3"

#Use BEDtools to identify the gene promoter (defined as 1000 bp upstream
#of a TSS) with the highest median signal in 'ctcf.hela.chr22.bg.gz'.
#Report the gene name (e.g., 'ABC123')

TSS="$datasets/bed/tss.hg19.chr22.bed.gz"

CTCF="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

hg19="$datasets/genome/hg19.genome"

answer_4=$(slopBed -i $TSS -g $hg19 -l 1000 -r 0 -s \
    | sortBed \
    | mapBed -a - -b $CTCF -c 4 -o median \
    | sort -k7nr \
    | cut -f 4 \
    | head -n 1)

echo "answer-4: $answer_4"

#Use BEDtools to identify the longest interval on 'chr22' that is not
#covered by 'genes.hg19.bed.gz'. Report the interval like 'chr1:100-500'.

genes="$datasets/bed/genes.hg19.bed.gz"

hg19="$datasets/genome/hg19.genome"

answer_5=$(sortBed -i $genes \
    | complementBed -i - -g $hg19 \
    | awk 'BEGIN {OFS="\t"} ($1=="chr22") {print $0, $3-$2}' \
    | sort -k4nr \
    | awk '(NR==1) {print $1":"$2"-"$3}')

echo "answer-5: $answer_5"

#(Extra Credit) Use one or more BEDtools that we haven't covered in class.

# EZH2 and SUZ12 are component of the Polycomb Repressor Complex 2, which has beed
# known to bind near or on CpG islands. I have decided to calculate the
# relative distance between overlapping sites where these two TFs bind and
# CpG islands using bedtools a

TFBS="$datasets/bed/encode.tfbs.chr22.bed.gz"

CpG="$datasets/bed/cpg.bed.gz"

gzcat $TFBS \
    | awk '($4 == "EZH2")' > ezh2.tfbs.bed

gzcat $TFBS \
    | awk '($4 == "SUZ12")' > suz12.tfbs.bed

answer_6=$(intersectBed -a ezh2.tfbs.bed -b suz12.tfbs.bed \
    | bedtools reldist -a - -b $CpG \
    | head -n 20)

rm ezh2.tfbs.bed suz12.tfbs.bed

echo -e "answer-6: \n$answer_6"
echo -e "\nAn over representation about 0.0 in the reldist column indicates
association of Ezh2/Suz12 complex and CpG islands."

