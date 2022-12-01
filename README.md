# 2022-sourmash-scaled-1

Calculates Jaccard similarity between two or more genomes, at full
resolution (no FracMinHash downsampling, equiv. scaled == 1).

Usage:
```
./jaccard-scaled-1.py genome1.fa genome2.fa.gz genome3.fq.gz -o numpy.mat

sourmash plot numpy.mat
```
