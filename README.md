## Summary
This repository contains two main files, which one should be able to execute as scripts and obtain all of the information asked for in the assigment.

### Part one (summarize partitions of a genome assembly)
The file "**partone.sh**" contains the complete script for the first part of the assignment. I've used "big" and "small" to refer to partitions with >1000 kb sequences and <1000 kb sequences respectively. 

Total number of nucleotides--
 * Larger sequence partition: 137547967
 * Smaller sequence partition: 6179905
 
Total number of Ns--
 * Larger sequence partition: 490385
 * Smaller sequence partition: 662593

Total number of sequences--
 * Larger sequence partition: 7
 * Smaller sequence partition: 1863

Nine of the .png files relate to this section. They are: 
1. GCall.png
2. GCbig.png
3. GCsmall.png
4. AllChr_CDF.png
5. big_CDF.png
6. small_CDF.png
7. lengthall.png
8. lengthbig.png
9. lengthsmall.png

__On a personal note__, I want to apologize for the quality of the histograms. For some reason, I couldn't get the R scripts that I was running from the command line to work with changes from defaults in the hist() function (axis labels, axis limits, and breaks), even though they were working in R studio on my personal computer.

### Part two (genome assembly)
The file "**parttwo.sh**" contains the script for the second part of the assignment.

For the assembly assessment section:
1. N50 of our nanopore assembly was 4494246, while the reference is 21485538
2. The dotplot comparing our assembly to the contig assembly from FlyBase is in the file named "dotplot.png"
3. The contiguity plot comparing the FlyBase assemblies with our assembly is in the file named "nanoporecomparison.png"
4. The BUSCO score for our assembly was 13/2799 (0.5%), whereas that of the reference was 2763/2799 (98.7%)

All of these files can be found in two different directories within your public directory: /pub/jje/ee282/galentm/HW4/ and /pub/jje/ee282/galentm/nanopore_assembly/

I haven't included the R scripts referenced by "partone.sh" because the repository already contains a lot of files, and I don't want to make it any harder to look through. They are located in /pub/jje/ee282/galentm/HW4/code/scripts/ 
Here is an example of one of them (they're all basically the same):
```
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
#sink(file = "/pub/jje/ee282/galentm/HW4/data/processed/GCall.png")
setwd("/pub/jje/ee282/galentm/HW4/data/processed/")

GCalldf <- read.table("allGC.txt", header = FALSE)
GCallvector <- GCalldf[[1]]

GCall.png <- hist(GCallvector, xlim = c(0,0.8), xlab = "GC content")
png("GCall.png")
plot(GCall.png)
dev.off()
```
