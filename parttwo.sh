#!/usr/bin/env bash
module load perl
module load jje/jjeutils
module load rstudio/0.99.9.9
module load jje/kent

#assembly code (aped entirely from you)
minimap=$(which minimap)
miniasm=$(which miniasm)
basedir=/pub/jje/ee282/$USER
projname=nanopore_assembly
basedir=$basedir/$projname
raw=$basedir/$projname/data/raw
processed=$basedir/$projname/data/processed
figures=$basedir/$projname/output/figures
reports=$basedir/$projname/output/reports

createProject $projname $basedir
ln -sf /bio/share/solarese/hw4/rawdata/iso1_onp_a2_1kb.fastq $raw/reads.fq

$minimap -t 32 -Sw5 -L100 -m0 $raw/reads.fq{,} \
| gzip -1 \
> $processed/onp.paf.gz

$miniasm -f $raw/reads.fq $processed/onp.paf.gz \
> $processed/reads.gfa

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| tee >(n50 /dev/stdin > $reports/n50.txt) \
| fold -w 60 \
> $processed/unitigs.fa

#split scaffold assembly into contigs
cd /pub/jje/ee282/galentm/nanopore_assembly/nanopore_assembly/data/raw/
faSplitByN dmel-all-chromosome-r6.24.fasta splitassembly.fa 10
mv splitassembly.fa /pub/jje/ee282/galentm/nanopore_assembly/nanopore_assembly/data/processed/

#N50 code (creates the function)
n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
}

#This chunk of code calculates N50s and stores them as text files. "flybase_n50.txt" can be compared to the value on NCBI to make sure it worked
cd ../processed/
n50 unitigs.fa > myassembly_n50.txt
n50 splitassembly.fa > flybase_n50.txt

#dotplot construction
###Loading of binaries via module load or PATH reassignment
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot/4.6.0

###Query and Reference Assignment. State my prefix for output filenames
REF="splitassembly.fa"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls u*.fa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)

###please use a value between 75-150 for -c. The value of 1000 is too strict.
nucmer -l 100 -c 75 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta \
  -R ${REF} -Q ${QRY} --png

#comparison between different assemblies with CDF
bioawk -c fastx ' { print length($seq) } ' splitassembly.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nFB_contig\t0" } { print "FB_contig\t" $1 } ' \
> FB_contig.tmp

bioawk -c fastx ' { print length($seq) } ' unitigs.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nanopore_contig\t0" } { print "nanopore_contig\t" $1 } ' \
> nanopore_contig.tmp

cd ../raw/
cp dmel-all-chromosome-r6.24.fasta /pub/jje/ee282/galentm/nanopore_assembly/nanopore_assembly/data/processed/
cd ../processed/

bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nFB_scaff\t0" } { print "FB_scaff\t" $1 } ' \
> FB_scaff.tmp

plotCDF2 {FB_scaff,FB_contig,nanopore_contig}.tmp nanoporecomparison.png

#BUSCO analysis
module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc
cd /pub/jje/ee282/galentm/nanopore_assembly/nanopore_assembly/data/processed/

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
##OPTIONS="${OPTIONS} -sp 4577"
QRY="unitigs.fa"
MYEXT=".fa" ###Please change this based on your qry file. I.e. .fasta or .fa or .gfa

#my busco run
BUSCO.py -c ${NSLOTS} -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}

module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc
cd /pub/jje/ee282/galentm/nanopore_assembly/nanopore_assembly/data/processed/

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
##OPTIONS="${OPTIONS} -sp 4577"
QRY="splitassembly.fa"
MYEXT=".fa" ###Please change this based on your qry file. I.e. .fasta or .fa or .gfa

#my busco run
BUSCO.py -c ${NSLOTS} -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}

