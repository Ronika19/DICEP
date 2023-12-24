#!/usr/bin/env perl 

use File::Copy;
use Statistics::R;
#use Bio::SeqIO;
use List::Util 'max';
use List::Util 'min';
use List::MoreUtils qw(uniq);
use experimental 'smartmatch';
use Getopt::Long qw(GetOptions);

OPTIONS();

# If output filename is provided by user
if (defined $outfile) {
	$output=shift(@ARGV); print("$output\n");
	$name=$output; 
}

#If sequence file is not provided
if (!defined $fasta) {
	if (scalar(@ARGV)<1){
	print("Sequence file not provided\n");
	exit;
	}
}

#If sequence file is provided
else {
	$f = shift @ARGV;
	$in_filename = $f;
	$name = $f if(!defined $outfile);
	$name =~ s/\..*//; #print("$name\n");
	#Check if sequence file exists
	unless (-e $f) {
		print("Sequence file does not exist\n");
		exit;
	} 
	$nce = 1;
	#Check if sequence file is in fasta format	
	open(F,$f);
	$line = <F>;
	if ($line !~ /^>/){
		print("Incorrect sequence file. Please provide sequence file in fasta format\n");
		exit;
	}
	#Get sequence from file
	while($line = <F>){
		$line =~ s/\n//;
		if($line !~ /^>/){
			$dna .= uc($line);
		}
	}
}

#Check only ATGC bases are present in the genome
@bases = ('A','T','G','C');
@genome = split("",$dna);
foreach(@genome){
	$s = $bases[rand@bases];
	$_ =~ s/[^ATGC]/$s/g;
}
if (!defined $verb) { $verb = 0; }
open (OUT, ">$name\_DICEP_1"); print("$name\n");
$m = join("",@genome);
print(OUT "$m");
$lenm = length($m);
print("Bacterial genome loaded\nLength of the genome is $lenm bp\n") if ($verb == 1); 

#If annotation file is provided
if (!defined $annotation){
	if (scalar(@ARGV) < 1){
		print("Annotation file does not exist\n");
		exit; 
	}

	$f1 = shift @ARGV;

	#Check if annotation file exists
	unless (-e $f1) {
		print("Annotation file does not exist\n");
		exit;
	} 
	
	open(F1,$f1);
	@annot_ptt = <F1>;
	
	#Check if annotation file is in ptt format
	$a = 0; $b = 0; $c = 0;
	foreach $annot_ptt(@annot_ptt) {
		@C = split(/\t/,$annot_ptt);
		$coords = $C[0];
		if (($coords =~ m/\.\./) && ($#C == 8)) {
			($start,$end) = split('\.\.',$coords); #print("$start\t$end\n");
			$arr1[$a] = $start;
			$arr2[$b] = $end;
			$arr3[$c] = $C[8]; # Gene Product
			$arr4[$d] = $C[3]; # PID
			$a++; $b++; $c++; $d++;
		}
	}
	
	@identifiers = split('\t',$annot_ptt[2]);

	if (defined $identifiers[0] && defined $identifiers[1] && defined $identifiers[8]) {
		if (($identifiers[0] !~ /Location/) | ($identifiers[1] !~ /Strand/) | ($identifiers[8] !~ /Product/)) {
			print("Incorrect annotation file. Please provide annotation file in ptt format\n");
			exit;
		}
	}

	else {
		print("Incorrect annotation file. Please provide annotation file in ptt format\n");
		exit;
	}

	#Get annotations from ptt file
	open (PTTOUT, ">$name\_gene_coord");
	$gene_counter = 0;
	for($i = 0; $i < $a; $i++){
		chomp $arr3[$i];
		print(PTTOUT "$arr3[$i]\t$arr1[$i]\t$arr2[$i]\n");
		$st = $arr1[$i]; $en = $arr2[$i]; $fn = $arr3[$i]; $prot_id = $arr4[$i];
		$gene_start[$gene_counter] = $st; $gene_end[$gene_counter] = $en; $gene_func[$gene_counter] = $fn;
		$pid[$gene_counter] = $prot_id;
		$gene_counter += 1;
	}
}

#If annotation file is not provided
else{
	#Run Prodigal for predicting genes
	print("Predicting genes using prodigal...\n") if ($verb == 1);
	#print "$in_filename\t$name\n";
	system("prodigal -a $name\.faa -f sco -i $in_filename -o $name\_sco_file.txt -q"); #bartonella =.fna file
	
	$sco = "$name\_sco_file.txt";
	open (S, $sco);
	$a = 0;$b = 0;$c = 0;
	while (<S>) {
		chomp $_;
		next if ($_ =~ m/#/);
		@coord = split('_', $_);
		$start[$a] = $coord[1];
		$end[$b] = $coord[2];
		$strand[$c] = $coord[3];
		$a++; $b++; $c++;
	}
	print("Identified $a genes\n") if ($verb == 1);
	print("Scanning genome for presence of genomic island specific marker genes...\n") if ($verb == 1);
	
	#Run HMMER for annnotating genes
	system("hmmscan -o $name\_stdout_out.txt --tblout $name\_table.txt CAFE_DB.hmm $name\.faa");

	$f = "$name\_table.txt";
	open(F,$f);


	#open hmmer output file
	while (<F>){
		chomp $_;
		next if $_ =~ /#/;
		@line = split('\s+',$_);	
		push(@id1,$line[2]);
		push (@evalue,$line[4]);
		$annot = join(' ', @line[18..$#line]);
		push(@func1, $annot);
	}
	
	#Get all markers with BLAST evalue < 0.01
	%hash = ();
	for ($i = 0; $i < scalar(@id1); $i++){
		next if ($id1[$i] eq $id1[$i-1]);
		@idd = split('_',$id1[$i]);
		if ($evalue[$i] < 0.01){  #Count only if evalue is less than 0.01.
			$id = $idd[1];
			$genome = $idd[0];
			$func = $func1[$i];

			$hash{$id} = $func unless (!defined $idd[1]);
		}
	
		#print("$id\t$func[$i]\n");
	}

	#Make ptt file for annotated genes
	open (OUT, ">$name.ptt");
	print(OUT "$genome\n");
	print(OUT "$c proteins\n");
	print(OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n");

	$marker_counter = 0;
	for ($i = 0; $i < $a; $i++){
		$num = $i+1;
		$length = int(($end[$i]-$start[$i])/3);
	
		if (exists($hash{$num})){
			print(OUT "$start[$i]..$end[$i]\t$strand[$i]\t$length\t-\t-\t-\t-\t-\t$hash{$num}\ transposase\n");
			$arr1[$i] = $start[$i]; $arr2[$i] = $end[$i]; $arr3[$i] = $hash{$num};
			$marker_counter += 1;
		}
		else {
			print(OUT "$start[$i]..$end[$i]\t$strand[$i]\t$length\t-\t-\t-\t-\t-\t-\n");
			$arr1[$i] = $start[$i]; $arr2[$i] = $end[$i]; $arr3[$i] = "-";
		}	
	}
	
	#open (PTTOUT, ">$name\_gene_coord");
	$gene_counter = 0;

	for($i = 0; $i < $a; $i++){
	#print(PTTOUT "$arr3[$i]\t$arr1[$i]\t$arr2[$i]\n");
	$st = $arr1[$i]; $en = $arr2[$i]; $fn = $arr3[$i];
	$gene_start[$gene_counter] = $st; $gene_end[$gene_counter] = $en; $gene_func[$gene_counter] = $fn;
	$gene_counter += 1;
	}
	print("Identified $marker_counter genomic island marker genes\n") if ($verb == 1);
}

#If phylo module is used
if ($phylo == 1){
	print("Checking phylogenetic distribution of genes\n") if ($verb == 1);
	print("$input_genus\n"); print("\n\n");
	#shift(@ARGV);
	system("cat ./faa/* > faa_database");
	@faa_files = <./faa/*>;
	$expect = scalar(@faa_files); #print("$expect\n");
	if ($expect < 1) {
		print("Phylo module requires atleast 1 genome for comparison\n");
		exit;
	}

	if(defined $annotation) {$infile = $name;}
	##make database
	system("makeblastdb -dbtype prot -in faa_database");#faa_database
	system("blastp -db faa_database -query $infile\.faa -outfmt \"6 qseqid sseqid stittle salltitles qcovs pident\" -out $name\_blast_output -num_threads 4") if(defined $annotation); 
	system("blastp -db faa_database -query $name\.faa -outfmt \"6 qseqid sseqid stittle salltitles qcovs pident\" -out $name\_blast_output -num_threads 4") if(!defined $annotation);  
	$blastf = "$name\_blast_output";
	open(BL, $blastf);
	open (BLOUT, ">$name\_phy"); open(PHYOUT, ">$name\_phyout");
	
	$bf = 0;
	while (<BL>){
		chomp;
		@BFC = split('\t', $_);
		$queryacc[$bf] = $BFC[0];
		#$subacc[$bf] = $BFC[1];
		$blsub[$bf] = $BFC[2];
		$coverage[$bf] = $BFC[3];
		$identities[$bf] = $BFC[4];
		$bf++;
	}
	$phprev="";
	$phcount=0;
	$phcum_count=0;
	$pi=0;
	foreach(@queryacc) {
		if($phprev ne $_) {
			if($phcount) {
				#printf("%s:%d \n",$_,$phcount);
				$phash{$phprev} = $phcount;
				$phash_cum{$phprev} = $phcum_count;
				#print("$_ $phcount $cum_count\n");
			}
			$ptemp_count = $phcount; #save values for last query
			$ptemp_cum_count = $phcum_count;
			$phprev = $_;
			$phcount = 0;
		}
		$phcount++;
		$phcum_count++;
		$pi++;
	}

	for ($i=0; $i<$bf; $i++) {
		if ($blsub[$i] =~ m/\[(.*)\]/){ #get blast subject names
			$blastname[$i]=$1; #print("$blastname[$i]\n");		
		}
	}

	for ($i=0; $i<$bf; $i++){
		$cur = $blastname[$i];
		if (defined $cur){
			@blast = split('\s+',$cur);
		}
		else { $blastname[$i] = ''; }
		$blast_genus[$i] = $blast[0];
		#$blast_sp[$i] = $blast[1];
		$blast_strain[$i] = join('', @blast[2..$#blast]);
		#print("$i $blast_strain[$i]\n");
		#print("$blast_genus[$i]\t$blast_sp[$i]\t$blast_strain[$i]\n");
	}
	
	@uniq_queryacc = uniq(@queryacc);

	# Extrating Protein Ids from faa file
	if (!defined $annotation) {
		for ($i=0; $i<scalar(@uniq_queryacc); $i++) {
			$protein_identifier[$i] = (split('\|', $uniq_queryacc[$i]))[1]; #print("$uniq_queryacc[$i]\t$protein_identifier[$i]\n");
		}
	}
	
	$last_element = $uniq_queryacc[-1];
	$phash{$last_element} = $ptemp_count;
	$phash_cum{$last_element} = $ptemp_cum_count;

	$spcount = 0;
	@all_strains = '';
	$phv = 0;
	
	for ($k=0; $k<scalar(@uniq_queryacc); $k++){
		$cur_query = $uniq_queryacc[$k];
		$next_query = $uniq_queryacc[$k+1]; #print("$cur_query\t$next_query\n");
		pop(@all_strains);
		$start = $phash_cum{$cur_query}-$phash{$cur_query};
		$end = $phash_cum{$cur_query};
		#print("$start\t$end\t$input_genus\n");
		for ($i=$start; $i<$end; $i++){
			
			if ($blast_genus[$i] =~ m/$input_genus/i){
				print(PHYOUT "$i\t$cur_query\t$phash{$cur_query}\t$coverage[$i]\t$identities[$i]\t$blast_strain[$i]\n");
				if (($coverage[$i] >= 70) && ($identities[$i] >= 70)) { #if identities and coverage are greater than 70
					$spcount += 1;
					push(@all_strains,$blast_strain[$i]);
					#print("$i\t$cur_query\t$spcount\t$phash{$cur_query}\t$blast_strain[$i]\n");
				}
			}
		}
		#print("@all_strains\n");
		$actual_count = scalar(uniq(@all_strains)); #Do not count multiple times if query matches multiple genes in same genome/strain
		
		if (!defined $expect ) {$pc=10000;}
		else { $pc = $actual_count/$expect; }
		if ($pc >= 0.5){
			$value = 0; #typical phyletic pattern
			$phyvalue[$phv] = 0;
			$phv++;
		}
		else {
			$value = 1; # atypical phyletic pattern
			$phyvalue[$phv] = 1;
			$phv++;
		} 
		
		if (!defined $expect) {
			print(BLOUT "$cur_query\t0\t$actual_count\t$phash{$cur_query}\t$value\n");
		}
		
		else {
			print(BLOUT "$cur_query\t$expect\t$actual_count\t$phash{$cur_query}\t$value\n");
		}
		
		$spcount=0;
		@all_strains='';
	}
}

#Run segmentation and clustering algorithm using compiled c program.
$segm=0.9999999999; $clust1=0.9999999999999; $clust2=0.9999999999999;

#LABEL:
#If user has defined thresholds
if (defined ($seg) && defined ($clus1) && defined ($clus2)){
	system("./segmentation_clustering.out $name\_DICEP_1 $seg $clus1 $clus2 $verb"); 
}

#Default thresholds
else {
	system("./segmentation_clustering.out $name\_DICEP_1 $segm $clust1 $clust2 $verb");	
} 	

#Parse segment and clustering output file 
$sc_out='segmentation_clustering.txt'; #seg-clus output file
open (SCOUT, $sc_out);
@SC=<SCOUT>;
open(F1,">$name\_DICEP_out");
$sca=0;$scb=0;$scc=0;
foreach (@SC){
	chomp;
	@SC_C=split(/\t/,$_); #print("$SC_C[0]\t$SC_C[1]\t$SC_C[2]\n");
	$sc_arr1[$sca++]=$SC_C[0];
	$sc_arr2[$scb++]=$SC_C[1];
	$sc_arr3[$scc++]=$SC_C[2];
}

#Combine contiguous segments with same cluster id
for ($i=0; $i<scalar(@sc_arr3); $i++) {
	if ($sc_arr3[$i] == $sc_arr3[$i+1]) {
		$sc_arr2[$i] = $sc_arr2[$i+1];
		$sc_arr1[$i+1] = $sc_arr1[$i]; 
	}
}

@unique_sc_arr3 = uniq(@sc_arr3);
@cond_arr1=(); @cond_arr2=(); @cond_arr3=(); @sc_arr=(); $x=0;
for ($i=0; $i<scalar(@sc_arr3); $i++) {
	#print(F1 "$sc_arr1[$i]\t$sc_arr2[$i]\t$sc_arr3[$i]\n");
	@indices = grep { $sc_arr1[$_] == $sc_arr1[$i] } 0..$#sc_arr1; # Find all the indices of duplicated elements in an array
	@duplicate_sc_arr2=(); #print("@indices\n");

	for ($j=0; $j<scalar(@indices); $j++) {
		push(@duplicate_sc_arr2, $sc_arr2[$indices[$j]]); 
	}

	if (! grep(/^$sc_arr1[$i]$/, @sc_arr)) {
		push(@sc_arr, $sc_arr1[$i]);
		push(@cond_arr1, $sc_arr1[$i]); 
		if (scalar(@indices) > 1) { push(@cond_arr2, max(@duplicate_sc_arr2)); }
		elsif (scalar(@indices) == 1) { push(@cond_arr2, $sc_arr2[$i]);}
		push(@cond_arr3, $sc_arr3[$i]);
		print(F1 "$cond_arr1[$x]\t$cond_arr2[$x]\t$cond_arr3[$x]\n");
		$x++;
	}
}

#Calculate size of each cluster as a percentage of the length of genome
open(F3,">$name\_DICEP_clustersize");
@unique_cond_arr3 = uniq(@cond_arr3); @perseg_endt = ();
for ($i=0; $i<scalar(@unique_cond_arr3); $i++) {
	$cluster_length[$i] = 0; 

	for ($j=0; $j<scalar(@cond_arr3); $j++) {
		if ($cond_arr3[$j] == $unique_cond_arr3[$i]) {
			$cluster_length[$i] += ($cond_arr2[$j]-$cond_arr1[$j]+1); 
		}
	}

	$perseg_endt[$i] = ($cluster_length[$i]/$lenm)*100;
	#print(F3 "$unique_cond_arr3[$i]\t$cluster_length[$i]\t$perseg_endt[$i]\n");
}

# Sorting the clusters in ascending order based on cluster number
%sort_hash=(); 
for ($i=min(@unique_cond_arr3); $i<=max(@unique_cond_arr3); $i++) {
	for ($j=0; $j<scalar(@unique_cond_arr3); $j++) {
		if ($unique_cond_arr3[$j] == $i) {
			print(F3 "$unique_cond_arr3[$j]\t$cluster_length[$j]\t$perseg_endt[$j]\n");
			$sort_hash{$cluster_length[$j]} = $unique_cond_arr3[$j];
		}
	}
}

# Sorting the clusters in ascending order based on cluster size
@sortedpc = sort { $a <=> $b } @perseg_endt;
$large_pc1 = $sortedpc[-1]; # Size percent of Largest cluster 
$large_pc2 = $sortedpc[-2]; # Size percent of 2nd Largest cluster
#print("$large_pc1\t$large_pc2\n"); 

@sorted_cluslen = sort {$a <=> $b} @cluster_length;
$large_cl1 = $sort_hash{$sorted_cluslen[-1]}; # Largest cluster
$large_cl2 = $sort_hash{$sorted_cluslen[-2]}; # 2nd Largest cluster
#print("$large_cl1\t$large_cl2\n");

######################## Using Composition Bias & Marker Enrichment OR Composition Bias & Aberrant Phyletic Pattern ######################

#Assign genes to clusters identified by segmentation and clustering algorithm
open (PYO1, ">$name\_PhyGenes");
$gex_counter=0;
for ($i=0; $i<scalar(@cond_arr3); $i++) {
	for ($j=0; $j<$gene_counter; $j++) {
		if (($cond_arr1[$i]<=$gene_start[$j]) && ($cond_arr1[$i]<$gene_end[$j]) && ($cond_arr2[$i]>$gene_start[$j]) && ($cond_arr2[$i]>=$gene_end[$j])) {
			$sc_clus[$gex_counter] = $cond_arr3[$i]; 
			$sc_start[$gex_counter] = $cond_arr1[$i]; $sc_end[$gex_counter] = $cond_arr2[$i]; 
			$g_start[$gex_counter] = $gene_start[$j]; $g_end[$gex_counter] = $gene_end[$j]; 
			$g_func[$gex_counter] = $gene_func[$j]; $g_func[$gex_counter] =~ s/[^a-zA-Z0-9]/ /g; 

			if (!defined $annotation) {
				if (grep(/^$pid[$j]$/, @protein_identifier)) {
					for ($x=0; $x<scalar(@protein_identifier); $x++) {
						if ($pid[$j] == $protein_identifier[$x]) {
							$phylovalue[$gex_counter] = $phyvalue[$x];
							print(PYO1 "$cond_arr3[$i]\t$cond_arr1[$i]\t$cond_arr2[$i]\t$gene_func[$j]\t$gene_start[$j]\t$gene_end[$j]\t$phyvalue[$x]\n");
							$gex_counter+=1;
						}	
					}
				}
			}
			
			elsif (defined $annotation) {
				$phylovalue[$gex_counter] = $phyvalue[$j]; 
				print(PYO1 "$cond_arr3[$i]\t$cond_arr1[$i]\t$cond_arr2[$i]\t$gene_func[$j]\t$gene_start[$j]\t$gene_end[$j]\t$phyvalue[$j]\n");		
				$gex_counter+=1;
			}
		}

		elsif (($cond_arr1[$i]<=$gene_start[$j]) && ($cond_arr1[$i]<$gene_end[$j]) && ($cond_arr2[$i]>$gene_start[$j]) && ($cond_arr2[$i]<=$gene_end[$j])) {
			$sc_clus[$gex_counter] = $cond_arr3[$i]; 
			$sc_start[$gex_counter] = $cond_arr1[$i]; $sc_end[$gex_counter] = $cond_arr2[$i]; 
			$g_start[$gex_counter] = $gene_start[$j]; $g_end[$gex_counter] = $gene_end[$j]; 
			$g_func[$gex_counter] = $gene_func[$j]; $g_func[$gex_counter] =~ s/[^a-zA-Z0-9]/ /g; 

			if (!defined $annotation) {
				if (grep(/^$pid[$j]$/, @protein_identifier)) {
					for ($x=0; $x<scalar(@protein_identifier); $x++) {
						if ($pid[$j] == $protein_identifier[$x]) {
							$phylovalue[$gex_counter] = $phyvalue[$x];
							print(PYO1 "$cond_arr3[$i]\t$cond_arr1[$i]\t$cond_arr2[$i]\t$gene_func[$j]\t$gene_start[$j]\t$gene_end[$j]\t$phyvalue[$x]\n");
							$gex_counter+=1;
						}	
					}
				}
			}
			
			elsif (defined $annotation) {
				$phylovalue[$gex_counter] = $phyvalue[$j]; 
				print(PYO1 "$cond_arr3[$i]\t$cond_arr1[$i]\t$cond_arr2[$i]\t$gene_func[$j]\t$gene_start[$j]\t$gene_end[$j]\t$phyvalue[$j]\n");		
				$gex_counter+=1;
			}
		}

		elsif (($cond_arr1[$i]>=$gene_start[$j]) && ($cond_arr1[$i]<$gene_end[$j]) && ($cond_arr2[$i]>$gene_start[$j]) && ($cond_arr2[$i]>=$gene_end[$j])) {
			$sc_clus[$gex_counter] = $cond_arr3[$i]; 
			$sc_start[$gex_counter] = $cond_arr1[$i]; $sc_end[$gex_counter] = $cond_arr2[$i]; 
			$g_start[$gex_counter] = $gene_start[$j]; $g_end[$gex_counter] = $gene_end[$j]; 
			$g_func[$gex_counter] = $gene_func[$j]; $g_func[$gex_counter] =~ s/[^a-zA-Z0-9]/ /g; 

			if (!defined $annotation) {
				if (grep(/^$pid[$j]$/, @protein_identifier)) {
					for ($x=0; $x<scalar(@protein_identifier); $x++) {
						if ($pid[$j] == $protein_identifier[$x]) {
							$phylovalue[$gex_counter] = $phyvalue[$x];
							print(PYO1 "$cond_arr3[$i]\t$cond_arr1[$i]\t$cond_arr2[$i]\t$gene_func[$j]\t$gene_start[$j]\t$gene_end[$j]\t$phyvalue[$x]\n");
							$gex_counter+=1;
						}	
					}
				}
			}
			
			elsif (defined $annotation) {
				$phylovalue[$gex_counter] = $phyvalue[$j]; 
				print(PYO1 "$cond_arr3[$i]\t$cond_arr1[$i]\t$cond_arr2[$i]\t$gene_func[$j]\t$gene_start[$j]\t$gene_end[$j]\t$phyvalue[$j]\n");		
				$gex_counter+=1;
			}
		}
	}
}

# Identifying Marker genes in each cluster
@Dlib=qw(transposase transposon integrase integration phage prophage bacteriophage mobile mobility insertion recombinase plasmid resolvase);
open(PYOUT2, ">$name\_ME_Genes");
@clus_gcounter = (); @total_counter = (); @marker_genes_start = (); @marker_genes_end = (); @gcounter = (); $total_marker_gene_counter = 0; 
for ($i=0; $i<scalar(@unique_cond_arr3); $i++) { # For each cluster
	for ($j=0; $j<scalar(@Dlib); $j++) { # For each marker term in @Dlib
		for ($k=0; $k<$gex_counter; $k++) { # For each gene
			$counter = 0; 
			if ($sc_clus[$k] == $unique_cond_arr3[$i]) {
				@split_g_func = split('\s+', $g_func[$k]); #print(@split_g_func);
				foreach (@split_g_func) { # For each word in gene function
					if (($_ ~~ /$Dlib[$j]/i) && (length($_) == length($Dlib[$j]))) {
						$counter += 1; # Marker term counter in a gene
						
						if ($gcounter[$k] < 1) {
							$gcounter[$k] += 1; 
						}
						if ((! grep(/^$g_start[$k]$/, @marker_genes_start)) && (! grep(/^$g_end[$k]$/, @marker_genes_end))) {
							push(@marker_genes_start, $g_start[$k]); 
							push(@marker_genes_end, $g_end[$k]);
							print(PYOUT2 "$unique_cond_arr3[$i]\t$g_start[$k]\t$g_end[$k]\t$counter\n");
						}
					}
				}
			}
			$total_counter[$k] = $counter; # Storing Total Marker terms in each gene
			if ($gcounter[$k] != 1) { $gcounter[$k] = 0; }
			@split_g_func = (); 
		}
	}
	foreach (@gcounter) {
		# Marker gene count in each cluster
		if ($_ == 1) { 
			$clus_gcounter[$unique_cond_arr3[$i]] += 1; 
			$total_marker_gene_counter += 1; # Total no of Marker genes in all clusters 
		}
		else { $clus_gcounter[$unique_cond_arr3[$i]] += 0; }
	}
	#print("$unique_cond_arr3[$i]\t$clus_gcounter[$unique_cond_arr3[$i]]\n");
	@gcounter = ();
}


$non_marker_genes_clusters = $gene_counter-$total_marker_gene_counter; # Genes that are not Markers (Non-Marker genes) in all the clusters
# Hypergeometric test for each cluster with Marker gene(s)
for ($i=0; $i<scalar(@unique_cond_arr3); $i++) { # For each cluster
	$cluster_gene_counter = 0;
	for ($k=0; $k<$gex_counter; $k++) { # For each gene
		if ($sc_clus[$k] == $unique_cond_arr3[$i]) {
			$cluster_gene_counter += 1;
		}
	}
	$non_marker_cluster[$unique_cond_arr3[$i]] = $cluster_gene_counter-$clus_gcounter[$unique_cond_arr3[$i]];
	$clus_genecounter[$unique_cond_arr3[$i]] = $cluster_gene_counter; # Total no of genes in a cluster set1
	$noncluster_gene[$unique_cond_arr3[$i]] = $gex_counter-$clus_genecounter[$unique_cond_arr3[$i]]; # Total no of genes not in a cluster, i.e., total no of non-cluster genes
		
	if (($clus_genecounter[$unique_cond_arr3[$i]] > 0) and ($clus_gcounter[$unique_cond_arr3[$i]] > 0)){
		$htest = hypergeometric($clus_gcounter[$unique_cond_arr3[$i]], $total_marker_gene_counter, $noncluster_gene[$unique_cond_arr3[$i]], $clus_genecounter[$unique_cond_arr3[$i]]);
		if ($htest > 0.5) { $sig[$unique_cond_arr3[$i]] = 1-$htest; }
		elsif ($htest <= 0.5) { $sig[$unique_cond_arr3[$i]] = $htest; }
		#print("Marker Info: $unique_cond_arr3[$i]\t$clus_gcounter[$unique_cond_arr3[$i]]\t$sig[$unique_cond_arr3[$i]]\n");
	}
	elsif (($clus_genecounter[$unique_cond_arr3[$i]] == 0) or ($clus_gcounter[$unique_cond_arr3[$i]] == 0)){
		$sig[$unique_cond_arr3[$i]] = 1;
	}
}


# Aberrant Phyletic Pattern & Marker Enrichment analysis
open(MEABB, ">$name\_ME_ABB");
$clus_gene_count = 0; $phylo_count = 0;
for ($i=0; $i<scalar(@unique_cond_arr3); $i++) {
	for ($j=0; $j<$gex_counter; $j++) {

		# Counting the number of genes in each cluster
		if ($sc_clus[$j] == $unique_cond_arr3[$i]) {
			$clus_gene_count += 1;

			# Counting the number of genes with aberrant phyletic pattern
			if ($phylovalue[$j] == 1) {
				$phylo_count += 1;
			}
			else { $phylo_count += 0; }
		}
	}

	# Aberrant Phyletic pattern analysis of each cluster. You need to choose a theshold for aberrant phyletic pattern analysis.
	if ($phylo_count > 0) {
		$aberrant = $phylo_count/$clus_gene_count;
		if ($aberrant >= 0.9) { # Aberrant Phyletic pattern threshold 
			$aberrant_cluster[$unique_cond_arr3[$i]] = $unique_cond_arr3[$i]; # Storing cluster with aberrant phyletic pattern
			$aberrant_phylo_value[$unique_cond_arr3[$i]] = $aberrant; # Storing aberrant phyletic values
		}
	}
	else { $aberrant = 0; }

	# Marker Enrichment analysis of each cluster. You need to choose a theshold for Marker Enrichment.
	if ($clus_gcounter[$unique_cond_arr3[$i]] > 0) {
		$marker_percentage = $clus_gcounter[$unique_cond_arr3[$i]]/$clus_gene_count;
		#if ($marker_percentage > 0.085) { # Marker enrichment threshold 
		$observed_markers_cluster = $clus_gcounter[$unique_cond_arr3[$i]]; 
		$expected_markers_cluster = $total_marker_gene_counter*($clus_gene_count/$gene_counter); 
		$fold_change_markers_cluster[$unique_cond_arr3[$i]] = $observed_markers_cluster/$expected_markers_cluster; 
		if (($sig[$unique_cond_arr3[$i]] <= 0.05) && ($fold_change_markers_cluster[$unique_cond_arr3[$i]] >= 4.5)) { 
			#if (($sig[$unique_cond_arr3[$i]] <= 0.05)) {
			$enrich[$unique_cond_arr3[$i]] = $marker_percentage; # Storing marker enrichment values
			$enriched_cluster[$unique_cond_arr3[$i]] = $unique_cond_arr3[$i]; # Storing marker enriched clusters
			#}
		}
	}
	else { $marker_percentage = 0; }
	print(MEABB "$unique_cond_arr3[$i]\t$clus_gene_count\t$marker_percentage\t$phylo_count\t$aberrant\t$fold_change_markers_cluster[$unique_cond_arr3[$i]]\t$sig[$unique_cond_arr3[$i]]\n"); 
	$clus_gene_count = 0; $phylo_count = 0;
}

# Determine Native & Alien clusters
push(@large_cl, $large_cl1);
push(@large_cl, $large_cl2);

@all_cond_arr3 = @cond_arr3;

open(NatAl, ">$name\_Native_Alien"); open(NatAlM, ">$name\_Native_Alien_Marker"); open(NatAlP, ">$name\_Native_Alien_Phylogenetic");
$x = 0; $y = 0;
for ($j=0; $j<scalar(@unique_cond_arr3); $j++) {
	for ($i=0; $i<scalar(@all_cond_arr3); $i++) {
		if ($all_cond_arr3[$i] == $unique_cond_arr3[$j]) {

			if (grep(/^$all_cond_arr3[$i]$/, @large_cl)) {
				if (! grep(/^$all_cond_arr3[$i]$/, @enriched_cluster)) {
					$native_cond_arr3[$x] = $all_cond_arr3[$i];
					$native_cond_arr1[$x] = $cond_arr1[$i]; 
					$native_cond_arr2[$x] = $cond_arr2[$i];
					$all_cond_arr3[$i] = 111111; # Native
					print(NatAl "$unique_cond_arr3[$j]\t$native_cond_arr3[$x]\t$native_cond_arr1[$x]\t$native_cond_arr2[$x]\t$all_cond_arr3[$i]\tNative\n");
					$x += 1;
				}
				elsif (grep(/^$all_cond_arr3[$i]$/, @enriched_cluster)) {
					$alien_cond_arr3[$y] = $all_cond_arr3[$i];
					$alien_cond_arr1[$y] = $cond_arr1[$i]; 
					$alien_cond_arr2[$y] = $cond_arr2[$i];
					$all_cond_arr3[$i] = 222222; # Alien
					print(NatAl "$unique_cond_arr3[$j]\t$alien_cond_arr3[$y]\t$alien_cond_arr1[$y]\t$alien_cond_arr2[$y]\t$all_cond_arr3[$i]\tAlien\n");
					$y += 1;
				}
			}

			elsif (! grep(/^$all_cond_arr3[$i]$/, @large_cl)) {

				if (grep(/^$all_cond_arr3[$i]$/, @enriched_cluster)) {
					if (grep(/^$all_cond_arr3[$i]$/, @aberrant_cluster)) {
						$alien_cond_arr3[$y] = $all_cond_arr3[$i];
						$alien_cond_arr1[$y] = $cond_arr1[$i]; 
						$alien_cond_arr2[$y] = $cond_arr2[$i];
						$all_cond_arr3[$i] = 222222; # Alien
						print(NatAl "$unique_cond_arr3[$j]\t$alien_cond_arr3[$y]\t$alien_cond_arr1[$y]\t$alien_cond_arr2[$y]\t$all_cond_arr3[$i]\tAlien\t$sig[$unique_cond_arr3[$j]]\n");
						print(NatAlM "$unique_cond_arr3[$j]\t$alien_cond_arr3[$y]\t$alien_cond_arr1[$y]\t$alien_cond_arr2[$y]\t$all_cond_arr3[$i]\tAlien\t$sig[$unique_cond_arr3[$j]]\n"); 
						$y += 1;
					}
					else {
						$alien_cond_arr3[$y] = $all_cond_arr3[$i];
						$alien_cond_arr1[$y] = $cond_arr1[$i]; 
						$alien_cond_arr2[$y] = $cond_arr2[$i];
						$all_cond_arr3[$i] = 222222; # Alien
						print(NatAl "$unique_cond_arr3[$j]\t$alien_cond_arr3[$y]\t$alien_cond_arr1[$y]\t$alien_cond_arr2[$y]\t$all_cond_arr3[$i]\tAlien\t$sig[$unique_cond_arr3[$j]]\n");
						print(NatAlM "$unique_cond_arr3[$j]\t$alien_cond_arr3[$y]\t$alien_cond_arr1[$y]\t$alien_cond_arr2[$y]\t$all_cond_arr3[$i]\tAlien\t$sig[$unique_cond_arr3[$j]]\n"); 
						$y += 1;
					}
				}

				elsif (! grep(/^$all_cond_arr3[$i]$/, @enriched_cluster)) {
					if (grep(/^$all_cond_arr3[$i]$/, @aberrant_cluster)) {
						$alien_cond_arr3[$y] = $all_cond_arr3[$i];
						$alien_cond_arr1[$y] = $cond_arr1[$i]; 
						$alien_cond_arr2[$y] = $cond_arr2[$i];
						$all_cond_arr3[$i] = 222222; # Alien
						print(NatAl "$unique_cond_arr3[$j]\t$alien_cond_arr3[$y]\t$alien_cond_arr1[$y]\t$alien_cond_arr2[$y]\t$all_cond_arr3[$i]\tNative\n");
						print(NatAlP "$unique_cond_arr3[$j]\t$alien_cond_arr3[$y]\t$alien_cond_arr1[$y]\t$alien_cond_arr2[$y]\t$all_cond_arr3[$i]\tNative\n"); 
						$y += 1;
					}
					else {
						$native_cond_arr3[$x] = $all_cond_arr3[$i];
						$native_cond_arr1[$x] = $cond_arr1[$i]; 
						$native_cond_arr2[$x] = $cond_arr2[$i];
						$all_cond_arr3[$i] = 111111; # Native
						print(NatAl "$unique_cond_arr3[$j]\t$native_cond_arr3[$x]\t$native_cond_arr1[$x]\t$native_cond_arr2[$x]\t$all_cond_arr3[$i]\tNative\n");
						$x += 1;
					}
				}
			}
		}
	}
}

# Sorting the Alien Clusters in ascending order based on cluster start & end
@sorted_alien_cond_arr1 = sort { $a <=> $b} @alien_cond_arr1;
@sorted_alien_cond_arr2 = sort { $a <=> $b} @alien_cond_arr2;

# Merging Contiguous Alien clusters
for ($i=0; $i<scalar(@sorted_alien_cond_arr1); $i++) {
	if ($i > 0) {
		if ($sorted_alien_cond_arr1[$i]-$sorted_alien_cond_arr2[$i-1] == 1) {
			$sorted_alien_cond_arr1[$i] = $sorted_alien_cond_arr1[$i-1]; 
			#print("$sorted_alien_cond_arr1[$i]\t$sorted_alien_cond_arr2[$i]\n");
		}
	}
}

$j = 0; @alien_cluster_start = (); @alien_cluster_end = ();
for ($i=0; $i<scalar(@sorted_alien_cond_arr1); $i++) {
	# Find all the indices of duplicated elements in an array
	@duplicate_index = grep { $sorted_alien_cond_arr1[$_] == $sorted_alien_cond_arr1[$i] } 0..$#sorted_alien_cond_arr1; 
	foreach (@duplicate_index) {
		push(@all_alien_cond_arr2, $sorted_alien_cond_arr2[$_]);
	}
	if (! grep(/^$sorted_alien_cond_arr1[$i]$/, @alien_cluster_start)) {
		# These are final Alien Clusters. These clusters needs to filtered based on their size i.e., >=8000 bp.
		$alien_cluster_start[$j] = $sorted_alien_cond_arr1[$i];
		$alien_cluster_end[$j] = max(@all_alien_cond_arr2); #print("$alien_cluster_start[$j]\t$alien_cluster_end[$j]\n"); 
		$j++;
	}
	@all_alien_cond_arr2 = ();
}

################################################## Using Aberrant Phyletic Pattern Only ##################################################

# Extract segments that display Aberrant Phyletic pattern
$y = 0;
for ($i=0; $i<$gene_counter; $i++) {
	$count_abb_genes = 0; $count_seg_genes = 0;
	if (defined $annotation) { 
		$phylogenetic_value = 0; 
	}
	elsif (!defined $annotation) { 
		if (grep(/^$pid[$i]$/, @protein_identifier)) {
			for ($x=0; $x<scalar(@protein_identifier); $x++) {
				if ($pid[$i] == $protein_identifier[$x]) {
					$phylogenetic_value = $phyvalue[$x];
				}
			}
		}
	}
	if ($phylogenetic_value == 1) {
		for ($j=0; $j<$gene_counter; $j++) {
			# Check to see if the gene has 4000 kb upstream and downstream
			if (($gene_start[$i] >= 4000) && (($gene_end[$i] >= 8000) or ($gene_end[$i] <= 8000))) {
				if ($lenm-$gene_end[$i] > 4000) {
					$gi_seg_start = $gene_start[$i]-4000;
					$gi_seg_end = $gene_end[$i]+4000;
					if (($gi_seg_start <= $gene_start[$j]) and ($gi_seg_start < $gene_end[$j]) and ($gi_seg_end > $gene_start[$j]) and ($gi_seg_end >= $gene_end[$j])) {
						if (defined $annotation) { $phy_value = $phyvalue[$j]; }
						elsif (!defined $annotation) { 
							if (grep(/^$pid[$j]$/, @protein_identifier)) {
								for ($x=0; $x<scalar(@protein_identifier); $x++) {
									if ($pid[$j] == $protein_identifier[$x]) {
										$phy_value = $phyvalue[$x];
									}
								}						
							}
						}
						$count_seg_genes += 1;
						if ($phy_value == 1) {
							$count_abb_genes += 1;
						}
					}
				}
				elsif ($lenm-$g_end[$i] < 4000) {
					$gi_seg_end = $lenm;
					$gi_seg_start = $gene_start[$i]-(8000-($gi_seg_end-$gene_end[$i]));
					if (($gi_seg_start <= $gene_start[$j]) and ($gi_seg_start < $gene_end[$j]) and ($gi_seg_end > $gene_start[$j]) and ($gi_seg_end >= $gene_end[$j])) {
						if (defined $annotation) { $phy_value = $phyvalue[$j]; }
						elsif (!defined $annotation) { 
							if (grep(/^$pid[$j]$/, @protein_identifier)) {
								for ($x=0; $x<scalar(@protein_identifier); $x++) {
									if ($pid[$j] == $protein_identifier[$x]) {
										$phy_value = $phyvalue[$x];
									}
								}						
							}
						}
						$count_seg_genes += 1;
						if ($phy_value == 1) {
							$count_abb_genes += 1;
						}
					}
				}
			}
			# If the gene has less than 4000 kb upstream and has 4000 kb downstream
			elsif (($gene_start[$i] <= 4000) && (($gene_end[$i] <= 4000) or ($gene_end[$i] >= 4000))) {
				$gi_seg_start = 1;
				$gi_seg_end = $gene_end[$i]+(8000-($gene_start[$i]-$gi_seg_start));
				if (($gi_seg_start <= $gene_start[$j]) and ($gi_seg_start < $gene_end[$j]) and ($gi_seg_end > $gene_start[$j]) and ($gi_seg_end >= $gene_end[$j])) {
					if (defined $annotation) { $phy_value = $phyvalue[$j]; }
					elsif (!defined $annotation) { 
						if (grep(/^$pid[$j]$/, @protein_identifier)) {
							for ($x=0; $x<scalar(@protein_identifier); $x++) {
								if ($pid[$j] == $protein_identifier[$x]) {
									$phy_value = $phyvalue[$x];
								}
							}						
						}
					}
					$count_seg_genes += 1;
					if ($phy_value == 1) {
						$count_abb_genes += 1;
					}
				}
			}
		}
		
		# Determining Thresholds for genes displaying Aberrant Phyletic pattern
		if (($count_abb_genes > 0) && ($count_seg_genes > 0)) {
			$abb_gene_percent = $count_abb_genes/$count_seg_genes;
			if ($abb_gene_percent > 0.99) { 
				# Storing the segment coordinates, gene counts and percent aberrant phyletic pattern displaying genes
				$abb_gene_count[$y] = $count_abb_genes; $seg_gene_count[$y] = $count_seg_genes;
				$gi_seg_st[$y] = $gi_seg_start; $gi_seg_en[$y] = $gi_seg_end; $abb_gi_percent[$y] = $abb_gene_percent;
				$y += 1; 
				#print("$gi_seg_start\t$gi_seg_end\t$count_abb_genes\t$count_seg_genes\t$abb_gene_percent\n");
			}
		}
		$count_abb_genes = 0; $count_seg_genes = 0;
	}
}

# Merging segments showing aberrant phyletic pattern
for ($i=0; $i<scalar(@gi_seg_st); $i++) {
	if ($i > 0) {
		if (($gi_seg_st[$i] < $gi_seg_en[$i-1]) && ($gi_seg_en[$i] > $gi_seg_en[$i-1])) {
			$gi_seg_st[$i] = $gi_seg_st[$i-1]; #print("$gi_seg_st[$i]\t$gi_seg_en[$i]\n");			
		}
	}
}

open(ABB_GI, ">$name\_ABB_GI.txt");
$p = 0; @seg_abb_start = (); @seg_abb_end = ();
for ($i=0; $i<scalar(@gi_seg_st); $i++) {
	@dup_seg_index = grep { $gi_seg_st[$_] == $gi_seg_st[$i] } 0..$#gi_seg_st; # Find all the indices of duplicated elements in an array
	foreach (@dup_seg_index) {
		push(@seg_abb_end, $gi_seg_en[$_]);
	}
	if (! grep(/^$gi_seg_st[$i]$/, @seg_abb_start)) {
		# These are final Alien segments displaying aberrant phyletic pattern. These segments needs to be filtered based on their size i.e., >=8000 bp.
		$segment_abb_start[$p] = $gi_seg_st[$i];
		$segment_abb_end[$p] = max(@seg_abb_end); 
		$p++;
	}
	@seg_abb_end = ();
}	

@uniqe_segment_abb_start = uniq(@segment_abb_start); @uniqe_segment_abb_end = uniq(@segment_abb_end);

# Writing Aberrant phyletic pattern displaying GIs in a file
$k = 0;
for ($i=0; $i<scalar(@uniqe_segment_abb_start); $i++) {
	if ($uniqe_segment_abb_end[$i]-$uniqe_segment_abb_start[$i] >= 8000) { 
		$final_abb_segment_start[$k] = $uniqe_segment_abb_start[$i];
		$final_abb_segment_end[$k] = $uniqe_segment_abb_end[$i];
		$segment_abb_length = $uniqe_segment_abb_end[$i]-$uniqe_segment_abb_start[$i];
		print(ABB_GI "$uniqe_segment_abb_start[$i]\t$uniqe_segment_abb_end[$i]\t$segment_abb_length\n"); 
		$k++;
	}
}

###################################################### Using Marker Enrichment Only ######################################################

# Extract segments that display Marker Enrichment
$y = 0;
for ($i=0; $i<$gene_counter; $i++) {
	$count_marker_genes = 0; $count_segment_genes = 0; $me_counter = 0;
	@split_gene_func = split('\s+', $gene_func[$i]); #print("$gene_func[$i]\n");
	for ($m=0; $m<scalar(@Dlib); $m++) {
		foreach (@split_gene_func) { # For each word in gene function
			if (($_ ~~ /$Dlib[$m]/i) && (length($_) == length($Dlib[$m]))) {
				$me_counter += 1; #print("$me_counter\n"); # Marker term counter in a gene
			}
		}
	}
	#print("$gene_func[$i]\t$me_counter\n");

	if ($me_counter > 0) {
		for ($j=0; $j<$gene_counter; $j++) {
			# Check to see if the gene has 4000 kb upstream and downstream
			if (($gene_start[$i] >= 4000) && (($gene_end[$i] >= 8000) or ($gene_end[$i] <= 8000))) {
				if ($lenm-$gene_end[$i] > 4000) {
					$gi_me_seg_start = $gene_start[$i]-4000;
					$gi_me_seg_end = $gene_end[$i]+4000;
					if (($gi_me_seg_start <= $gene_start[$j]) and ($gi_me_seg_start < $gene_end[$j]) and ($gi_me_seg_end > $gene_start[$j]) and ($gi_me_seg_end >= $gene_end[$j])) {
						@split_me_gene_func = split('\s+', $gene_func[$j]); 
						$counter = 0;
						for ($m=0; $m<scalar(@Dlib); $m++) {
							foreach (@split_me_gene_func) { # For each word in gene function
								if (($_ ~~ /$Dlib[$m]/i) && (length($_) == length($Dlib[$m]))) {
									$counter += 1; # Marker term counter in a gene 
									#print("$gi_me_seg_start\t$gi_me_seg_end\t$gene_start[$j]\t$gene_end[$j]\t$gene_func[$j]\t$counter\n");
								}
							}
						}
						if ($counter > 0) { $counter = 1;  }
						elsif ($counter == 0) { $counter = 0; }
						$count_segment_genes += 1;
						$count_marker_genes += $counter; 
						#print("$gi_me_seg_start\t$gi_me_seg_end\t$gene_start[$j]\t$gene_end[$j]\t$gene_func[$j]\t$count_segment_genes\t$count_marker_genes\n");
					}
				}

				elsif ($lenm-$g_end[$i] < 4000) {
					$gi_me_seg_end = $lenm;
					$gi_me_seg_start = $gene_start[$i]-(8000-($gi_me_seg_end-$gene_end[$i]));
					if (($gi_me_seg_start <= $gene_start[$j]) and ($gi_me_seg_start < $gene_end[$j]) and ($gi_me_seg_end > $gene_start[$j]) and ($gi_me_seg_end >= $gene_end[$j])) {
						@split_me_gene_func = split('\s+', $gene_func[$j]); 
						$counter = 0;
						for ($m=0; $m<scalar(@Dlib); $m++) {
							foreach (@split_me_gene_func) { # For each word in gene function
								if (($_ ~~ /$Dlib[$m]/i) && (length($_) == length($Dlib[$m]))) {
									$counter += 1; # Marker term counter in a gene 
									#print("$gi_me_seg_start\t$gi_me_seg_end\t$gene_start[$j]\t$gene_end[$j]\t$gene_func[$j]\t$counter\n");
								}
							}
						}
						if ($counter > 0) { $counter = 1;  }
						elsif ($counter == 0) { $counter = 0; }
						$count_segment_genes += 1;
						$count_marker_genes += $counter; 
						#print("$gi_me_seg_start\t$gi_me_seg_end\t$gene_start[$j]\t$gene_end[$j]\t$gene_func[$j]\t$count_segment_genes\t$count_marker_genes\n");
					}
				}
			}

			# If the gene has less than 4000 kb upstream and has 4000 kb downstream
			elsif (($gene_start[$i] <= 4000) && (($gene_end[$i] <= 4000) or ($gene_end[$i] >= 4000))) {
				$gi_me_seg_start = 1;
				$gi_me_seg_end = $gene_end[$i]+(8000-($gene_start[$i]-$gi_me_seg_start));
				if (($gi_me_seg_start <= $gene_start[$j]) and ($gi_me_seg_start < $gene_end[$j]) and ($gi_me_seg_end > $gene_start[$j]) and ($gi_me_seg_end >= $gene_end[$j])) {
					@split_me_gene_func = split('\s+', $gene_func[$j]); 
					$counter = 0;
					for ($m=0; $m<scalar(@Dlib); $m++) {
						foreach (@split_me_gene_func) { # For each word in gene function
							if (($_ ~~ /$Dlib[$m]/i) && (length($_) == length($Dlib[$m]))) {
								$counter += 1; # Marker term counter in a gene 
								#print("$gi_me_seg_start\t$gi_me_seg_end\t$gene_start[$j]\t$gene_end[$j]\t$gene_func[$j]\t$counter\n");
							}
						}
					}
					if ($counter > 0) { $counter = 1;  }
					elsif ($counter == 0) { $counter = 0; }
					$count_segment_genes += 1;
					$count_marker_genes += $counter; 
					#print("$gi_me_seg_start\t$gi_me_seg_end\t$gene_start[$j]\t$gene_end[$j]\t$gene_func[$j]\t$count_segment_genes\t$count_marker_genes\n");
				}
			}
		}

		# Determining Thresholds for genes displaying Marker Enrichment
		if (($count_marker_genes > 0) && ($count_segment_genes > 0)) {
			$marker_gene_percent = $count_marker_genes/$count_segment_genes; 
			#print("$gi_me_seg_start\t$gi_me_seg_end\t$count_marker_genes\t$count_segment_genes\t$marker_gene_percent\n");
			$non_marker_count = $count_segment_genes-$count_marker_genes; # Non-Marker gene count in segment 
			$non_segment_genes = $gene_counter-$count_segment_genes; # No of genes not in a segment i.e., non-segmental genes 
			$htest_segment = hypergeometric($count_marker_genes, $total_marker_gene_counter, $non_segment_genes, $count_segment_genes); 
			if ($htest_segment > 0.5) { $sigval = 1-$htest_segment; } 
			elsif ($htest_segment <= 0.5) { $sigval = $htest_segment; } 

			$observed_markers_segment = $count_marker_genes; 
			$expected_markers_segment = $total_marker_gene_counter*($count_segment_genes/$gene_counter); 
			if ($expected_markers_segment == 0) { $expected_markers_segment = 1; }
			$fold_change_markers_segment = $observed_markers_segment/$expected_markers_segment; 
			if (($sigval <= 0.05) && ($fold_change_markers_segment >= 16.5)) { 
				# Storing the segment coordinates, gene counts and marker gene percent
				$marker_gene_count[$y] = $count_marker_genes; $segment_gene_count[$y] = $count_segment_genes;
				$gi_meseg_st[$y] = $gi_me_seg_start; $gi_meseg_en[$y] = $gi_me_seg_end; 
				$marker_gi_percent[$y] = $marker_gene_percent; 
				#print("$gi_me_seg_start\t$gi_me_seg_end\t$count_marker_genes\t$count_segment_genes\t$marker_gene_percent\t$fold_change_markers_segment\t$sigval\n"); 
				$y += 1; 
			}
		}

		$count_segment_genes = 0; $count_marker_genes = 0;
	}
}

# Merging segments showing marker enrichment
for ($i=0; $i<scalar(@gi_meseg_st); $i++) {
	if ($i > 0) {
		if (($gi_meseg_st[$i] < $gi_meseg_en[$i-1]) && ($gi_meseg_en[$i] > $gi_meseg_en[$i-1])) {
			$gi_meseg_st[$i] = $gi_meseg_st[$i-1]; #print("$gi_meseg_st[$i]\t$gi_meseg_en[$i]\n");			
		}
	}
}

open(ME_GI, ">$name\_ME_GI.txt");
$p = 0; @seg_me_start = (); @seg_me_end = ();
for ($i=0; $i<scalar(@gi_meseg_st); $i++) {
	@dup_meseg_index = grep { $gi_meseg_st[$_] == $gi_meseg_st[$i] } 0..$#gi_meseg_st; # Find all the indices of duplicated elements in an array
	foreach (@dup_meseg_index) {
		push(@seg_me_end, $gi_meseg_en[$_]);
	}
	if (! grep(/^$gi_meseg_st[$i]$/, @seg_me_start)) {
		# These are final Alien segments displaying marker enrichment. These segments needs to be filtered based on their size i.e., >=8000 bp.
		$segment_me_start[$p] = $gi_meseg_st[$i];
		$segment_me_end[$p] = max(@seg_me_end); 
		$p++;
	}
	@seg_me_end = ();
}

@uniqe_segment_me_start = uniq(@segment_me_start); @uniqe_segment_me_end = uniq(@segment_me_end);

# Writing Marker enriched GIs in a file
$k = 0;
for ($i=0; $i<scalar(@uniqe_segment_me_start); $i++) {
	if ($uniqe_segment_me_end[$i]-$uniqe_segment_me_start[$i] >= 8000) { 
		$final_me_segment_start[$k] = $uniqe_segment_me_start[$i];
		$final_me_segment_end[$k] = $uniqe_segment_me_end[$i];
		$segment_me_length = $uniqe_segment_me_end[$i]-$uniqe_segment_me_start[$i];
		print(ME_GI "$uniqe_segment_me_start[$i]\t$uniqe_segment_me_end[$i]\t$segment_me_length\n"); 
		$k++;
	}
}

# Storing the segments identified to be enriched in only Marker genes & the segments that display only Aberrant Phyletic pattern.
push(@final_abb_me_segment_start, @final_abb_segment_start); push(@final_abb_me_segment_start, @final_me_segment_start); 
push(@final_abb_me_segment_end, @final_abb_segment_end); push(@final_abb_me_segment_end, @final_me_segment_end);

# Sorting the Segments that are enriched in Markers & segments displaying Aberrant Phyletic pattern.
@sorted_final_abb_me_segment_start = sort { $a <=> $b } @final_abb_me_segment_start;
@sorted_final_abb_me_segment_end = sort { $a <=> $b } @final_abb_me_segment_end;

for ($i=0; $i<scalar(@sorted_final_abb_me_segment_start); $i++) {
	if ($i>0) {
		if (($sorted_final_abb_me_segment_start[$i] < $sorted_final_abb_me_segment_end[$i-1]) && ($sorted_final_abb_me_segment_end[$i] >= $sorted_final_abb_me_segment_end[$i-1])) {
			$sorted_final_abb_me_segment_start[$i] = $sorted_final_abb_me_segment_start[$i-1]; 
			#print("$sorted_final_abb_me_segment_start[$i]\t$sorted_final_abb_me_segment_end[$i]\n");
		}
	}
}

$k = 0; @sort_final_abb_me_segment_start = (); @sort_final_abb_me_segment_end = (); print("\n");
for ($i=0; $i<scalar(@sorted_final_abb_me_segment_start); $i++) {
	# Find all the indices of duplicated elements in an array
	@duplicate_indexes = grep { $sorted_final_abb_me_segment_start[$_] == $sorted_final_abb_me_segment_start[$i] } 0..$#sorted_final_abb_me_segment_start; 
	foreach (@duplicate_indexes) {
		push(@all_final_abb_me_end, $sorted_final_abb_me_segment_end[$_]);
	}
	if (! grep(/^$sorted_final_abb_me_segment_start[$i]$/, @sort_final_abb_me_segment_start)) {
		# These are final Marker Enriched segments & final Aberrant Phyletic pattern displaying segments. These segments needs to filtered based on their size i.e., >=8000 bp.
		$sort_final_abb_me_segment_start[$k] = $sorted_final_abb_me_segment_start[$i];
		$sort_final_abb_me_segment_end[$k] = max(@all_final_abb_me_end); 
		#print("$sort_final_abb_me_segment_start[$k]\t$sort_final_abb_me_segment_end[$k]\n"); 
		$k++;
	}
	@all_final_abb_me_end = ();
}

########################################## Using Aberrant Phyletic pattern & Marker Enrichment ###########################################

# Extract the segments from the Native clusters that display Aberrant Phyletic pattern & Marker Enrichment
@unique_native_cond_arr3 = uniq(@native_cond_arr3);
open(LC1 ,">$name\_LC1"); open(LC2, ">$name\_LC2");
$x = 0; @seg_start = (); @seg_end = ();
for ($k=0; $k<scalar(@unique_native_cond_arr3); $k++) {
	# Find if native clusters are enriched in markers or not. Just to double check that these clusters are native.
	if (! grep(/^$unique_native_cond_arr3[$k]$/, @enriched_cluster)) { # If the native clusters are not enriched in markers
		
		for ($j=0; $j<$gex_counter; $j++) {
			if ($sc_clus[$j] == $unique_native_cond_arr3[$k]) { # For these native clusters
				$count_genes = 0; $count_abb_me_genes = 0; $count_abb = 0; $count_me_genes = 0;
				for ($i=0; $i<$gex_counter; $i++) {
					if ($sc_clus[$i] == $unique_native_cond_arr3[$k]) { # For these largest clusters
						# Check to see if the gene has 4000 kb upstream and downstream
						if (($g_start[$j] >= 4000) && (($g_end[$j] >= 8000) or ($g_end[$j] <= 8000))) {
							if ($lenm-$g_end[$j] > 4000) {
								$gstart = $g_start[$j]-4000;
								$gend = $g_end[$j]+4000; 
								if ($gstart < $sc_start[$j]) {
									$gstart = $sc_start[$j];
								}
								if ($gend > $sc_end[$j]) {
									$gend = $sc_end[$j];
								}
								if (($gstart <= $g_start[$i]) and ($gstart < $g_end[$i]) and ($gend > $g_start[$i]) and ($gend >= $g_end[$i])) {
									print(LC1 "$sc_clus[$i]\t$gstart\t$gend\t$g_start[$i]\t$g_end[$i]\n");
									$count_genes += 1; 
	
									# Extract the segments in these clusters that display Aberrant Phyletic pattern & Marker Enrichment
									if ($phylovalue[$i] == 1) { 
										$count_abb += 1;
										if ((grep (/^$g_start[$i]$/, @marker_genes_start)) && (grep(/^$g_end[$i]$/, @marker_genes_end))) {	
											$count_abb_me_genes += 1; 
										}
									}
									if ((grep (/^$g_start[$i]$/, @marker_genes_start)) && (grep(/^$g_end[$i]$/, @marker_genes_end))) {
										$count_me_genes += 1;
									}
								}
							}

							elsif ($lenm-$g_end[$j] < 4000) {
								$gend = $lenm;
								$gstart = $g_start[$j]-(8000-($gend-$g_end[$j])); 
								if ($gstart < $sc_start[$j]) {
									$gstart = $sc_start[$j];
								}
								if ($gend > $sc_end[$j]) {
									$gend = $sc_end[$j];
								}
								if (($gstart <= $g_start[$i]) and ($gstart < $g_end[$i]) and ($gend > $g_start[$i]) and ($gend >= $g_end[$i])) {
									print(LC1 "$sc_clus[$i]\t$gstart\t$gend\t$g_start[$i]\t$g_end[$i]\n");
									$count_genes += 1; 

									# Extract the segments in these clusters that display Aberrant Phyletic pattern & Marker Enrichment
									if ($phylovalue[$i] == 1) {
										$count_abb += 1;
										if ((grep (/^$g_start[$i]$/, @marker_genes_start)) && (grep(/^$g_end[$i]$/, @marker_genes_end))) {	
											$count_abb_me_genes += 1; 
										}
									}
									if ((grep (/^$g_start[$i]$/, @marker_genes_start)) && (grep(/^$g_end[$i]$/, @marker_genes_end))) {
										$count_me_genes += 1;
									}
								}
							}
						}

						# If the gene has less than 4000 kb upstream and has 4000 kb downstream
						elsif (($g_start[$j] <= 4000) && (($g_end[$j] <= 4000) or ($g_end[$j] >= 4000))) {
							$gstart = 1;
							$gend = $g_end[$j]+(8000-($g_start[$j]-$gstart)); 
							if ($gstart < $sc_start[$j]) {
								$gstart = $sc_start[$j];
							}
							if ($gend > $sc_end[$j]) {
								$gend = $sc_end[$j];
							}
							if (($gstart <= $g_start[$i]) and ($gstart < $g_end[$i]) and ($gend > $g_start[$i]) and ($gend >= $g_end[$i])) {
								print(LC1 "$sc_clus[$i]\t$gstart\t$gend\t$g_start[$j]\t$g_end[$i]\n");
								$count_genes += 1; 

								# Extract the segments in these clusters that display Aberrant Phyletic pattern & Marker Enrichment
								if ($phylovalue[$i] == 1) {
									$count_abb += 1;
									if ((grep (/^$g_start[$i]$/, @marker_genes_start)) && (grep(/^$g_end[$i]$/, @marker_genes_end))) {												
										$count_abb_me_genes += 1; 
									}	
								}
								if ((grep (/^$g_start[$i]$/, @marker_genes_start)) && (grep(/^$g_end[$i]$/, @marker_genes_end))) {
									$count_me_genes += 1;
								}
							}
						}
					}
				}
				# Thresholds for Marker Enrichment and Aberrant Phyletic pattern
				if (($count_abb > 0) && ($count_me_genes > 0) && ($count_genes > 0)) {
					$me_genes_percentage = $count_me_genes/$count_genes;
					$abb_genes_percentage = $count_abb/$count_genes;
					if (($abb_genes_percentage >= 0.9) && ($me_genes_percentage > 0)) { 
						print(LC1 "$sc_clus[$j]\t$gstart\t$gend\t$count_genes\t$count_me_genes\t$count_abb\t$count_abb_me_genes\t$abb_genes_percentage\t$me_genes_percentage\n");
						if ((! grep(/^$gstart$/, @seg_start)) && (! grep(/^$gend$/, @seg_end))) {
							$seg_clus_id[$x] = $sc_clus[$j]; $seg_starts[$x] = $gstart; $seg_ends[$x] = $gend; 
							$seg_me_genes_percentage[$x] = $me_genes_percentage; $seg_abb_genes_percentage[$x] = $abb_genes_percentage;
							$seg_count_me_genes[$x] = $count_me_genes; 
							$seg_count_genes[$x] = $count_genes; 
							$nonseg_genes[$x] = $gene_counter-$count_genes; # No of genes not in the segment (Non-Segment genes)
							#print(LC2 "$clus_id[$x]\t$seg_start[$x]\t$seg_end[$x]\t$abb_genes_percentage\t$me_genes_percentage\n");
							$x++;
						}
					}
				}	
			}
		}
	}
}

for ($i=0; $i<scalar(@seg_starts); $i++) {
	$htest_segmt = hypergeometric($seg_count_me_genes[$i], $total_marker_gene_counter, $nonseg_genes[$i], $seg_count_genes[$i]); #print("$htest_segmt\n");
	if ($htest_segmt > 0.5) { $sigvalue = 1-$htest_segmt; } 
	elsif ($htest_segmt <= 0.5) { $sigvalue = $htest_segmt; }
	$observed_markers_segmt[$i] = $seg_count_me_genes[$i]; $expected_markers_segmt[$i] = $total_marker_gene_counter*($seg_count_genes[$i]/$gene_counter);
	if ($expected_markers_segment == 0) { $expected_markers_segment = 1; }
	$fold_change_markers_segmt[$i] = $observed_markers_segmt[$i]/$expected_markers_segmt[$i];
	if (($sigvalue <= 0.05) && ($fold_change_markers_segmt[$i] >= 4.5)) {
		$clus_id[$i] = $seg_clus_id[$i]; $seg_start[$i] = $seg_starts[$i]; $seg_end[$i] = $seg_ends[$i];
		print(LC2 "$clus_id[$i]\t$seg_start[$i]\t$seg_end[$i]\t$seg_abb_genes_percentage[$i]\t$seg_me_genes_percentage[$i]\t$fold_change_markers_segmt[$i]\t$sigvalue\n");
	}
}

# Merging Alien segments in Native Clusters
for ($i=0; $i<scalar(@seg_start); $i++) {
	if ($i > 0) {
		if (($seg_start[$i] < $seg_end[$i-1]) && ($seg_end[$i]) > $seg_end[$i-1]) {
			$seg_start[$i] = $seg_start[$i-1]; #print("$seg_start[$i]\t$seg_end[$i]\n");
		}
	}
}

$k = 0; @segment_start = (); @segment_end = ();
for ($i=0; $i<scalar(@seg_start); $i++) {
	@dup_index = grep { $seg_start[$_] == $seg_start[$i] } 0..$#seg_start; # Find all the indices of duplicated elements in an array
	foreach (@dup_index) {
		push(@all_segment_end, $seg_end[$_]);
	}
	if (! grep(/^$seg_start[$i]$/, @segment_start)) {
		# These are final Alien segments in Native Clusters. These segments needs to filtered based on their size i.e., >=8000 bp.
		$segment_start[$k] = $seg_start[$i];
		$segment_end[$k] = max(@all_segment_end); #print("$segment_start[$k]\t$segment_end[$k]\n"); 
		$k++;
	}
	@all_segment_end = ();
}				

# Getting all the Alien horizontally transferred Islands
@all_alien_start = (); @all_alien_end = ();
push(@all_alien_start, @segment_start); push(@all_alien_start, @alien_cluster_start);
push(@all_alien_end, @segment_end); push(@all_alien_end, @alien_cluster_end);

# Sorting the Alien Islands
@sort_all_alien_start = sort { $a <=> $b } @all_alien_start;
@sort_all_alien_end = sort { $a <=> $b } @all_alien_end;

# Determining and Writing all the Genomic Islands (Comp bias + ABB, Comp bias + ME, ABB + ME) to the output file
open(GI, ">$name\_GenomicIslands.txt");
for ($i=0; $i<scalar(@all_alien_start); $i++) {
	$alien_length = $sort_all_alien_end[$i]-$sort_all_alien_start[$i];
	if ($alien_length >= 8000) { # filtering Alien Islands based on size ####
		print(GI "$sort_all_alien_start[$i]\t$sort_all_alien_end[$i]\t$alien_length\n");
	}
}

##########################################################################################################################################

# Getting all the overlapping Genomic Islands ((Comp bias + ABB, Comp bias + ME, ABB + ME) & (ABB) & (ME))
$r = 0;
for ($i=0; $i<scalar(@all_alien_start); $i++) {
	$all_alien_length = $sort_all_alien_end[$i]-$sort_all_alien_start[$i];
	if ($all_alien_length >= 8000) { # filtering Alien Islands based on size 
		for ($j=0; $j<scalar(@sort_final_abb_me_segment_start); $j++) {
			if ($sort_final_abb_me_segment_end[$j]-$sort_final_abb_me_segment_start[$j] >= 8000) { 
				if (($sort_all_alien_start[$i] <= $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_start[$i] < $sort_final_abb_me_segment_end[$j]) && ($sort_all_alien_end[$i] > $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_end[$i] >= $sort_final_abb_me_segment_end[$j])) {
					$alien_gi_st[$r] = $sort_all_alien_start[$i]; 
					$alien_gi_en[$r] = $sort_all_alien_end[$i];
					$abb_gi_st[$r] = $sort_final_abb_me_segment_start[$j]; 
					$abb_gi_en[$r] = $sort_final_abb_me_segment_end[$j];
					$final_alien_start[$r] = $sort_all_alien_start[$i]; 
					$final_alien_end[$r] = $sort_all_alien_end[$i];
					#print("$sort_all_alien_start[$i]\t$sort_all_alien_end[$i]\t$sort_final_abb_me_segment_start[$j]\t$sort_final_abb_me_segment_end[$j]\t$final_alien_start[$r]\t$final_alien_end[$r]\n");
					$r++; 
				}
				elsif (($sort_all_alien_start[$i] <= $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_start[$i] < $sort_final_abb_me_segment_end[$j]) && ($sort_all_alien_end[$i] > $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_end[$i] <= $sort_final_abb_me_segment_end[$j])) {
					$alien_gi_st[$r] = $sort_all_alien_start[$i]; 
					$alien_gi_en[$r] = $sort_all_alien_end[$i];
					$abb_gi_st[$r] = $sort_final_abb_me_segment_start[$j]; 
					$abb_gi_en[$r] = $sort_final_abb_me_segment_end[$j];
					$final_alien_start[$r] = $sort_all_alien_start[$i]; 
					$final_alien_end[$r] = $sort_final_abb_me_segment_end[$j];
					#print("$sort_all_alien_start[$i]\t$sort_all_alien_end[$i]\t$sort_final_abb_me_segment_start[$j]\t$sort_final_abb_me_segment_end[$j]\t$final_alien_start[$r]\t$final_alien_end[$r]\n");
					$r++; 
				}
				elsif (($sort_all_alien_start[$i] >= $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_start[$i] < $sort_final_abb_me_segment_end[$j]) && ($sort_all_alien_end[$i] > $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_end[$i] >= $sort_final_abb_me_segment_end[$j])) {
					$alien_gi_st[$r] = $sort_all_alien_start[$i]; 
					$alien_gi_en[$r] = $sort_all_alien_end[$i];
					$abb_gi_st[$r] = $sort_final_abb_me_segment_start[$j]; 
					$abb_gi_en[$r] = $sort_final_abb_me_segment_end[$j];
					$final_alien_start[$r] = $sort_final_abb_me_segment_start[$j]; 
					$final_alien_end[$r] = $sort_all_alien_end[$i];
					#print("$sort_all_alien_start[$i]\t$sort_all_alien_end[$i]\t$sort_final_abb_me_segment_start[$j]\t$sort_final_abb_me_segment_end[$j]\t$final_alien_start[$r]\t$final_alien_end[$r]\n");
					$r++; 
				}
				elsif (($sort_all_alien_start[$i] >= $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_start[$i] < $sort_final_abb_me_segment_end[$j]) && ($sort_all_alien_end[$i] > $sort_final_abb_me_segment_start[$j]) && ($sort_all_alien_end[$i] <= $sort_final_abb_me_segment_end[$j])) {
					$alien_gi_st[$r] = $sort_all_alien_start[$i]; 
					$alien_gi_en[$r] = $sort_all_alien_end[$i];
					$abb_gi_st[$r] = $sort_final_abb_me_segment_start[$j]; 
					$abb_gi_en[$r] = $sort_final_abb_me_segment_end[$j];
					$final_alien_start[$r] = $sort_final_abb_me_segment_start[$j]; 
					$final_alien_end[$r] = $sort_final_abb_me_segment_end[$j];
					#print("$sort_all_alien_start[$i]\t$sort_all_alien_end[$i]\t$sort_final_abb_me_segment_start[$j]\t$sort_final_abb_me_segment_end[$j]\t$final_alien_start[$r]\t$final_alien_end[$r]\n");
					$r++; 
				}
			}
		}
	}
}

# Merging and Writing all the Genomic Islands (Comp bias + ABB, Comp bias + ME, ABB + ME, ABB, ME) to the output file
for ($i=0; $i<scalar(@final_alien_start); $i++) {
	$final_length[$i] = $final_alien_end[$i]-$final_alien_start[$i]; 
	#print("$final_alien_start[$i]\t$final_alien_end[$i]\n");
	if ($i > 0) {
		if (($final_alien_start[$i] < $final_alien_end[$i-1]) && ($final_alien_end[$i] >= $final_alien_end[$i-1])) {
			$final_alien_start[$i] = $final_alien_start[$i-1]; 
			#print("$final_alien_start[$i]\t$final_alien_end[$i]\n");
		}
	}
}

$k = 0; @final_island_start = (); @final_island_end = (); #print("\n");
for ($i=0; $i<scalar(@final_alien_start); $i++) {
	@dup_indexes = grep { $final_alien_start[$_] == $final_alien_start[$i] } 0..$#final_alien_start; # Find all the indices of duplicated elements in an array
	foreach (@dup_indexes) {
		push(@all_final_end, $final_alien_end[$_]);
	}
	if (! grep(/^$final_alien_start[$i]$/, @final_island_start)) {
		# These are final Alien segments and Alien Clusters. These segments needs to filtered based on their size i.e., >=8000 bp.
		$final_island_start[$k] = $final_alien_start[$i];
		$final_island_end[$k] = max(@all_final_end); 
		#print("$final_island_start[$k]\t$final_island_end[$k]\n"); 
		$k++;
	}
	@all_final_end = ();
}

# Removing duplicate genomic islands
$m = 0; $n = 0; #print("\n");
for ($i=0; $i<scalar(@sort_all_alien_start); $i++) {
	if (! grep(/^$sort_all_alien_start[$i]$/, @alien_gi_st)) {
		if ($sort_all_alien_end[$i]-$sort_all_alien_start[$i] >= 8000) { 
			$final_alien_island_start[$m] = $sort_all_alien_start[$i]; $final_alien_island_end[$m] = $sort_all_alien_end[$i];
			push(@final_island_start, $final_alien_island_start[$m]); push(@final_island_end, $final_alien_island_end[$m]);
			$m++; #print("$sort_all_alien_start[$i]\t$sort_all_alien_end[$i]\n");
		}
	}
}
print("\n");
for ($i=0; $i<scalar(@sort_final_abb_me_segment_start); $i++) {
	if (! grep(/^$sort_final_abb_me_segment_start[$i]$/, @abb_gi_st)) {
		if ($sort_final_abb_me_segment_end[$i]-$sort_final_abb_me_segment_start[$i] >= 8000) {
			$final_alien_island_start[$n] = $sort_final_abb_me_segment_start[$i]; $final_alien_island_end[$n] = $sort_final_abb_me_segment_end[$i];
			push(@final_island_start, $final_alien_island_start[$n]); push(@final_island_end, $final_alien_island_end[$n]);
			$n++;# print("$sort_final_abb_me_segment_start[$i]\t$sort_final_abb_me_segment_end[$i]\n");
		}
	}
}

@sorted_final_island_start = sort {$a <=> $b} @final_island_start;
@sorted_final_island_end = sort {$a <=> $b} @final_island_end;

# Writing genomic islands that were determined based on Comp+ME, Comp+ABB, ABB+ME, ABB, ME.
open(ABB_ME_COMP_ABB_ME, ">$name\_DICEP_All_GIs.txt");
for ($i=0; $i<scalar(@sorted_final_island_start); $i++) {
	$island_length[$i] = $sorted_final_island_end[$i]-$sorted_final_island_start[$i];
	print(ABB_ME_COMP_ABB_ME "$sorted_final_island_start[$i]\t$sorted_final_island_end[$i]\t$island_length[$i]\n");
}


################################################################## CGView to generate GI images ##########################################

# Generating genomic island map using CGview
if ($visual==1){
	print("Preparing files for making visualizing genomic islands\n") if ($verb==1);
	open(VOUT, ">$name\_DICEP_feature_table"); print(VOUT "seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n");

	for ($i=0; $i<scalar(@sorted_final_island_start); $i++){
		print VOUT "gene$i\t.\tgene\t$sorted_final_island_start[$i]\t$sorted_final_island_end[$i]\t1\t.\t.\n";
	}

	open(VISL, ">$name\_DICEP_label"); 
	print(VISL "gene\t1");
	if ($verb==0){
		system ("perl ~/cgview-master/scripts/cgview_xml_builder/cgview_xml_builder.pl -sequence $in_filename -verbose s -genes $name\_DICEP_feature_table -title \"$name\" -labels_to_show $name\_DICEP_label -custom featureThickness=45 labelFontSize=45 -output $name\_DICEP_marker.xml");
	}
	else {
		system ("perl ~/cgview-master/scripts/cgview_xml_builder/cgview_xml_builder.pl -sequence $in_filename -genes $name\_DICEP_feature_table -title \"$name\" -labels_to_show $name\_DICEP_label -custom featureThickness=45 labelFontSize=45 -output $name\_DICEP_marker.xml");
	}

	push (@vers, "$name\_DICEP_marker.xml");
	$ci=0;

	#Edit XML file to include labels and color GIs
	foreach $vers (@vers){ 
		open(VISF, $vers);
		@file=<VISF>;
		$count=0;
		open(VISF1,">$name\_DICEP1\_$ci.xml");	
		

		foreach(@file){
			if ($_=~m/showLabel/){
				$count+=1;
			}
	
		}		

		foreach(@file){
			if ($_=~m/showLabel/){
				$_=~s/showLabel="false"/label=\"GI-$count" showLabel="true"/;
				$count-=1;
			}
			if ($_=~m/plain, 20/){
				$_=~s/plain, 20/plain, 45/;
			}
			if ($_=~m/showShading="true"/){
				$_=~s/showShading="true"/showShading="false"/;
			}
			if ($_=~m/51,51,51/){
					$_=~s/51,51,51/0,0,204/;
			}
			if ($_=~m/plain, 80" text/){
					$_=~s/plain, 80" text/italics, 80" text/;
			}


			print VISF1 "$_";
		}				

		$ci+=1;
	}
	system ("java -jar ~/cgview-master/bin/cgview.jar -i $name\_DICEP1_0.xml -o $name\_DICEP_marker.png > $name\_DICEP_CGViewout");
	
	
}

######################################################################################################################################################

if ($expert==0){
	system ("rm $name\_stdout_out.txt") if (defined$annotation);
	system ("rm $name\_DICEP_1");	
	system ("rm $name\_DICEP_clustersize");
	system ("rm $name\_ABB_GI.txt");
	system ("rm $name\_DICEP_out");
	system ("rm $name\_gene_coord");
	system ("rm $name\_GenomicIslands.txt");
	system ("rm $name\_LC1");
	system ("rm $name\_LC2");
	system ("rm $name\_ME_ABB");
	system ("rm $name\_ME_Genes");
	system ("rm $name\_ME_GI.txt");
	system ("rm $name\_Native_Alien");
	system ("rm $name\_Native_Alien_Marker");
	system ("rm $name\_Native_Alien_Phylogenetic");
	system ("rm $name\_PhyGenes");
	system ("rm $name\_phyout");
	system ("rm segmentation_clustering.txt");
	system ("rm segments.txt");
	system ("rm contclus_round1.txt");
	system ("rm contclus_round2.txt");
	system ("rm $name\_phy") if($phylo==1); 
	system ("rm $name\_blast_output") if($phylo==1);	
	system ("rm faa_database") if (defined$phylo);	
	system ("rm $name\_sco_file.txt") if (defined$annotation);
	system ("rm $name\_table.txt") if (defined$annotation);
}


if (defined $annotation) {
	print("Completed scanning for genomic island specific marker genes!\n"); 
}

sub information {
	print("DICEP is a genomic island prediction tool. This is DICEP v1.0 2023, written by Ronika De <ronikade\@my.unt.edu>\n");
	exit;
}

sub usage {
	print(STDERR "\nUsage: When annotation file is available:\n  perl $0 [options] --phylogenetic --genus Escherichia --fasta NC_004431.fna NC_004431.ptt\n\n");
	print("\nUsage: When annotation file is not available:\n  perl $0 [options] --phylogenetic --genus Escherichia --annotation --fasta NC_004431.fna\n\n");
  	exit;
}

sub thresholds {
	$seg=shift(@ARGV);
	$clus1=shift(@ARGV);
	$clus2=shift(@ARGV);
	if (defined ($seg) && defined ($clus1) && defined ($clus2)){
		if ($seg=~m/^-?\d+\.?\d*$/ && $clus1=~m/^-?\d+\.?\d*$/ && $clus2=~m/^-?\d+\.?\d*$/){ #matches float type numbers

			if ($seg<0||$seg>1||$clus1<0||$clus1>1||$clus2<0||$clus2>1){ #check if thresholds are in range
				print("Thresholds should be in range 0 to 1\n");
				exit;
			}
			else {		
				print("Thresholds are $seg $clus1 $clus2\n");
			}
		}
		else{
			print("Input Segmentation, Contiguous-Clustering & Noncontiguous-Clustering thresholds\n");
			exit;
		}	
	}
	else {
		print("Input Segmentation, Contiguous-Clustering & Noncontiguous-Clustering thresholds\n");
		exit;
	}
}

sub verbose {
	$verb=1; print("Verbose Activated!\n");
}

sub expert {
	$expert=1; print("Keeping temporary files!\n");
}

sub visual {
	$visual=1;
}

sub phylogenetic {
	$phylo=1; print("Phylogenetic Module Activated!\n");
}

sub genus {
	$input_genus=shift(@ARGV);
	#print "Input genus $input_genus\n";
	return $input_genus;
}

sub hypergeometric {
	($x,$m,$n,$k) = @_; #print("$x,$m,$n,$k\n");
	$R_hyper = Statistics::R->new();
	$R_hyper->run( qq`x = phyper($x,$m,$n,$k, lower.tail=FALSE)` );
	$squares = $R_hyper->get('x'); #print("$squares\n");
	return $squares;
}

sub OPTIONS {
	GetOptions(
	"help" =>\&usage,
	"information" =>\&information,
	"annotation" =>\$annotation,
	"fasta" =>\$fasta,
	"thresholds" =>\&thresholds,
	"genbank" =>\&genbank,
	"phylogenetic" =>\&phylogenetic,
	"genus" =>\&genus,
	"out" =>\$outfile,
	"verbose" =>\&verbose,
	"expert" =>\&expert,
	"visual" =>\&visual,
	);
}




