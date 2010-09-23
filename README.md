jsnprank
========

#### SNP ranking algorithm ####

### Description ###
SNPrank[1] is an eigenvector centrality algorithm that ranks the importance of 
single nucleotide polymorphisms (SNPs) in a genetic association interaction 
network (GAIN) [2]. Each SNP is ranked according to its overall contribution 
to the phenotype, including its main effect and second- and higher-order 
gene-gene interactions. 

This software was created as a bioinformatics tool for usage by our research 
group, [In Silico](http://insilico.utulsa.edu), as well as other researchers 
and interested parties.  

### Dependencies and Usage ###
jsnprank is developed and tested on 64-bit Linux (Ubuntu), but should work on any 
platform supported by Java (OpenJDK java version 1.6.0 tested).

jsnprank requires the [jblas](http://jblas.org) linear algebra library 
for matrix computations (tested with jblas version 1.1.1).

To compile the java classes:

    javac -cp jblas-1.1.1.jar:. *.java

Once the Java source files have been compiled, you can run SNPrank from 
command-line:

    java -cp jblas-1.1.1.jar:. snprank matrix.txt 0.85 output.txt

### Contributors ###
See AUTHORS file.

### References ###
[1]N.A. Davis, J.E. Crowe, Jr., N.M. Pajewski, and B.A. McKinney. Surfing a 
genetic association interaction network to identify modulators of antibody 
response to smallpox vaccine. Genes and Immunity, 2010, 
doi: 10.1038/gene.2010.3. [open access](http://www.nature.com/gene/journal/vaop/ncurrent/full/gene201037a.html)

[2]B.A. McKinney, J.Guo, J.E. Crowe, Jr., and D. Tian. Capturing the spectrum of 
interaction effects in genetic association studies by simulated evaporative 
cooling network analysis. PLoS Genetics 2009, 5(3): e1000432. 
doi:10.1371/journal.pgen.1000432. [open access](http://www.plosgenetics.org/article/info:doi/10.1371/journal.pgen.1000432)
