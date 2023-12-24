################################## Readme File ##############################

Requirements:
1. Ubuntu OS
2. Perl version >= 5.22.0 
3. gcc version: 7.4.0


In-built C code usage:
The DICEP.pl uses a pre-compiled C code. 
The C code (segmentation_clustering.c) can be compiled using the command: $gcc segmentation_clustering.c -o segmentation_clustering.out -lm



Perl libraries required:
1. List::MoreUtils
# Installation command:
$perl -MCPAN -e shell
cpan> install List::MoreUtils

2. Statistics::R 
# Installation command: sudo  perl -MCPAN -e 'install Statistics::R'



Additional Tool Requirements:
1. NCBI-BLAST+ version: 2.6.0+	
# Installation command: $sudo apt install ncbi-blast+

2. Prodigal version: 2.6.3	
# Installation command:  $sudo apt-get install -y prodigal 

3. HMMER version: 3.1b2	
# Installation command:  $sudo apt-get install -y hmmer 



To run the perl script and generate output you need to run the following command when annotation files are provided:
$perl DICEP.pl --phylogenetic --genus Escherichia --fasta NC_004431.fna NC_004431.ptt

To run the perl script and generate output you need to run the following command when annotation files are not available:
$perl DICEP.pl --phylogenetic --genus Escherichia --annotation --fasta NC_004431.fna

Input Files collected from GenBank for bacterial genome, for example:
When annotation files are available:
1. 'NC_004431.fna'
2. 'NC_004431.faa'
3. 'NC_004431.ptt'

When annotation files are not available (i.e., for draft genomes):
1. 'NC_004431.fna'

Output Files:
The 'NC_004431_DICEP_All_GIs.txt' file is the final output file containing the genomic islands (GIs).



