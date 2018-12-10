#!/usr/bin/env bash
module load perl
module load jje/jjeutils/0.1a
module load rstudio/0.99.9.9
module load jje/kent

#Create directory for HW4, part 1
createProject HW4 /pub/jje/ee282/galentm/

#
#partition into sequences less than and greater than 100kb
cd /pub/jje/ee282/galentm/HW4/data/raw/
bioawk -c fastx ' length($seq) > 100000 { print $seq} ' dmel-all-chromosome-r6.24.fasta > /pub/jje/ee282/galentm/HW4/data/processed/biggerseq.fasta
bioawk -c fastx ' length($seq) < 100000 { print $seq} ' dmel-all-chromosome-r6.24.fasta > /pub/jje/ee282/galentm/HW4/data/processed/smallerseq.fasta

cd ../processed/

#calculate number of nucleotides, number of Ns, and number of sequences for smaller partition
cat smallerseq.fasta \
| tee >( \
	wc -c > smaller_bitcount.txt ) \
| tee >( \
	wc -l > smaller_seqcount.txt) \
| tr -cd N \
| wc -c > smaller_Ncount.txt

#calculate number of nucleotides, number of Ns, and number of sequences for bigger partition
cat biggerseq.fasta \
| tee >( \
	wc -c > bigger_bitcount.txt ) \
| tee >( \
	wc -l > bigger_seqcount.txt) \
| tr -cd N \
| wc -c > bigger_Ncount.txt


cd ../raw/

#create temporary files for plotCDF2 input
bioawk -c fastx ' length($seq) > 100000 { print length($seq) } ' dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nOver100kb\t0" } { print "Over100kb\t" $1 } ' > /pub/jje/ee282/galentm/HW4/data/processed/biggerlengths.tmp

bioawk -c fastx ' length($seq) < 100000 { print length($seq) } ' dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nUnder100kb\t0" } { print "Under100kb\t" $1 } ' > /pub/jje/ee282/galentm/HW4/data/processed/smallerlengths.tmp

bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk 'BEGIN { print "Assembly\tLength\nAllChr\t0" } { print "AllChr\t" $1 } ' > /pub/jje/ee282/galentm/HW4/data/processed/AllChrlengths.tmp

#plot with plotCDF
cd ../processed
plotCDF2 AllChrlengths.tmp AllChr_CDF.png
plotCDF2 smallerlengths.tmp small_CDF.png
plotCDF2 biggerlengths.tmp big_CDF.png

#create files with GC content of sequences in partitions and whole undivided data
cd ../raw
bioawk -c fastx ' { print gc($seq) } ' dmel-all-chromosome-r6.24.fasta > /pub/jje/ee282/galentm/HW4/data/processed/allGC.txt
bioawk -c fastx ' length($seq) > 100000 { print gc($seq) } ' dmel-all-chromosome-r6.24.fasta > /pub/jje/ee282/galentm/HW4/data/processed/bigGC.txt
bioawk -c fastx ' length($seq) < 100000 { print gc($seq) } ' dmel-all-chromosome-r6.24.fasta > /pub/jje/ee282/galentm/HW4/data/processed/smallGC.txt

#turns temporary files from plotCDF into appropriate format for R histogram script (removes 0)
cd ../processed
tail -n +3 smallerlengths.tmp > smallerlengths2.tmp
tail -n +3 AllChrlengths.tmp > AllChrlengths2.tmp
tail -n +3 biggerlengths.tmp > biggerlengths2.tmp

#These six R scripts create histograms showing GC and sequence length distribution for the different files.
cd ../../code/scripts/
./GChistsmall.R
./GChistbig.R
./GChistall.R
./alllengthfreq.R
./smalllengthfreq.R
./biglengthfreq.R

cd ../../data/processed/
rm *tmp
