/*
This script is part of a genomic island prediction program, CAFE
It first divides the genomes into several segments, following multiple rounds of binary segmentation
Next, it perform agglomerative hierarchical clustering in two steps- contiguous clustering and non-contiguous clustering
*/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h> 
#include<time.h>
#include<stddef.h>
#include <errno.h>

// Global constants
const int ITMAX = 100;
const double EPS = 3.0e-7;
const double FPMIN = 1.0e-30; 

//Declare functions
double entropy(double *freq_oligos, double);
int segmentation(int *hash,int parent_start, int parent_end);
int clustering(int final_s1_start[], int final_s1_end[], int strlength);
float gammp(float , float );
float gammln(float xx);

//Global variables
extern int errno;
int errnum;
int hash[15000000], m=2, nm,all_kmer, count_of_hash, sg1=1, sg2, final_s1_start[900000], final_s1_end[900000], s2j=0,s2_start[900000], s2_end[900000],l=0,k1, dof,mmk,fgs=0,verb;
double thres, thres1, thres3;

int main(int argc, char **argv) {
	//Get thresholds
	thres=atof(*(&argv[2]));//Default segmentation threshold=0.8
	thres1=atof(*(&argv[3]));//Default contiguous clustering threshold=0.99999
	thres3=atof(*(&argv[4]));//Default non-contiguous clusterin threshold=0.999
	//printf ("%lf %lf %lf\n", thres, thres1, thres3);
	verb=atoi(*(&argv[5]));

	int i, j, id,sum, length,k=0,exp,  num,n, count=0,   c=0, q=0, count1=0, ascend,strlength=0 ;
	char *genome=malloc(15000000*sizeof(char));
	char *dna=malloc(100*sizeof(char));
	double ent;

	//Check if input file provided
	if (argc==1) { fprintf(stderr,"Did you forget to include filename?\n"); }

	FILE *input, *noutput;
	input = fopen(argv[1], "r");

	if (input==NULL){
		errnum = errno;
		//fprintf(stderr, "Value of errno: %d\n", errno);
		//perror("Error printed by perror");
		fprintf(stderr, "Incorrect filename: %s\n", strerror( errnum ));
	}

	//open the input file 
	else {

		fgets (genome, 15000000, input);
		fclose(input);
		length=strlen(genome);

		//Convert DNA sequence to a numeric sequence and store in int array
		//m=markov model order, 2. This converts all trinucleotides to numbers 1-64. DNA is now represented by numbers 1-64, each representing a rinucleotide
		nm=m+1;
		all_kmer=pow(4.0,nm);
		dof=(pow(4.0,m))*3;
		mmk=(4*4*(pow(4,m))-1);

		for(i=0; i<length-nm; i++) {
			sum = 0;
			strncpy(dna, genome+i,nm);
			for(j=0;j<nm;j++){
				num=pow(10, n);
                		if (dna[j]=='A'){id=0;}
	        		else if (dna[j]=='T'){id=1;}
	        		else if (dna[j]=='C'){id=2;}
	        		else if (dna[j]=='G'){id=3;}
	        		exp=nm-j-1;
	        		sum = sum + (id * pow(4.0,exp));
			}
			strlength+=1;
			hash[i+1] = sum+1;
			count_of_hash+=1; 
		}
		strlength+=2;

	}

	clustering(final_s1_start, final_s1_end, strlength);
	free(genome); free(dna);
}

/*****************************************************************************************************************************************/

int clustering(int final_s1_start[], int final_s1_end[], int strlength) {
	//Initialize variables
	int count=0, count1=0, i, j, k=0, q=0, ascend;
	int p2length,  s,final_groups_start[5000], final_groups_end[5000],switc=0, np2length,hcounter=1, len=1,switcarray[5000]; 
	double enta, entb, entab, jsd, pi1, pi2, sig =0.3,p1length1, total_length, newp1l=0.0,limit=0.9999999999999;
	double a, b, c1, d, beta, neff, arg, sx, smax, sarg, smax1; 
	//double thres1=0.9999999999999;
	double *freqab=malloc((mmk+1)*sizeof(double));
	double *freqa=malloc((mmk+1)*sizeof(double));
	double *freqb=malloc((mmk+1)*sizeof(double));
	double *tempafreq=malloc((mmk+1)*sizeof(double));
	double *cluslen=malloc(10000*sizeof(double));
	double *hash1=malloc(100000*sizeof(double));
	double *total_cluster_length=malloc(1000*sizeof(double));

	if (verb==1) {printf("\nStarting segmentation...\nThis may take some time\n");}
	//printf("\nstrlength is %d\n", strlength);

	//Start segmentation
	if (verb==1) {printf ("Segment_start\tSegment_end\tDivergence\n");}
	segmentation(hash,1,strlength);

	for(i=1; i<sg1-1; i++) {
		for (j=i+1; j<sg1; j++) {
			if (final_s1_start[i] > final_s1_start[j]) {
				ascend =  final_s1_start[i];
				final_s1_start[i] = final_s1_start[j];
                		final_s1_start[j] = ascend;
			}

			if (final_s1_end[i] > final_s1_end[j]) {
				ascend =  final_s1_end[i];
				final_s1_end[i] = final_s1_end[j];
		        	final_s1_end[j] = ascend;
			}
		}
	}

	//Print all segments
	FILE *output1;    
	output1 = fopen("segments.txt", "w");
	for(i=1; i<sg1; i++) {
		if (verb==1) { fprintf(output1,"from main final_sement_one%d start:%d end:%d\n",i, final_s1_start [i], final_s1_end[i]); }
	}

	if (verb==1) {printf("Clustering starts...\n");}
	if (verb==1) {printf("It is almost done now\n");}
	FILE *output2 = fopen("contclus_round1.txt", "w"); FILE *output2_1 = fopen("contclus_round2.txt", "w");
	//FILE *output2_2 = fopen("contclus_round1_cluslen", "w"); FILE *output2_3 = fopen("contclus_round2_cluslen", "w");
	fgs=1; sg1=sg1-1;

	for (i=1; i<sg1; i++) { // for all segments
		if (switc==0) { p1length1 = final_s1_end [i]-final_s1_start[i]+1-3; }
		p2length = final_s1_end [i+1]-final_s1_start[i+1]-1;
		np2length=p2length-1;
		total_length=p1length1+p2length;

 		for(k=1; k<all_kmer+1; k++) {
			if (switc==0) {
				//Get frequencies for segment 1  
				for(s=final_s1_start [i]; s<final_s1_end[i]-2; s++) {
					if(k==hash[s]) { count=count+1; }
				}
				freqa[k]=count;
			}
			//Get frequencies for segment 2 
			for(s=final_s1_start [i+1]; s<final_s1_end[i+1]-2; s++) {
				if(k==hash[s]) { count1=count1+1; }
			}
			freqb[k]=count1;

			//Get frequencies for the combined segment
			freqab[k]=freqa[k]+freqb[k];
			count=0; count1=0;
                   	}
		//Calculate entropy and divergence between segments  
		enta = -((entropy(freqa,p1length1))/(log(2.0)));
		entb = -((entropy(freqb,p2length))/(log(2.0)));
		entab = -((entropy(freqab,total_length))/(log(2.0)));
		pi1 = (double)p1length1/total_length;
		pi2 = (double)1-pi1;
		jsd = entab-(pi1*enta)-(pi2*entb);

		//Markov model parameters
		if (m==2) { a=2.39; b=-7.66;c1=0.0029; d=0.841; }
		else if (m==0) { a=2.7784; b=-7.97084; c1=0.0; d=0.80; }				
		else if (m==1) { a=2.543; b=-4.77; c1=0.0; d=0.848; }

		//Do significance testing							
		beta=(c1*log(total_length)) + d;
		neff=(a*log(total_length)) + b;
		sarg=(log(2.0)*total_length*jsd);
		sx=gammp(dof/2,sarg);	
		smax=sx;
		//printf("%d %d %d %d %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %.11lf %lf %lf %lf %lf %.15lf\n",final_s1_start[i],final_s1_end [i],final_s1_start[i+1],final_s1_end[i+1],entab, enta, p1length1, pi1, pi1*enta, entb,p2length,pi2, pi2*entb,total_length, jsd,beta,neff,sarg,sx,smax);

		//Merge segments if divergence is not significant
		if (smax<limit) {
			switc+=1;
			if(switc==1) { //If a segment merges with only one segment
				//Get frequencies of merged segment-clusters
				p1length1=(double)(total_length-1)/2;
				newp1l=total_length-1;

				for(k=1; k<all_kmer+1; k++){
					freqa[k]=(double)((freqab[k])/2);
					tempafreq[k]=freqab[k];
				}
			}

			else{ //If segment merges with more than one segment
				//Get frequencies of merged segment-clusters
				p1length1=(double)(newp1l+np2length)/(switc+1);
				newp1l=newp1l+p2length-1;

				for(k=1;k<all_kmer+1;k++){
					freqa[k]=(double)(tempafreq[k]+freqb[k])/(switc+1);
					tempafreq[k]=tempafreq[k]+freqb[k];
				}
			}

			final_s1_start[i+1]=final_s1_start[i];

			if (i+1==sg1){//For last segment
				final_groups_start[fgs]=final_s1_start[i];
				final_groups_end[fgs]=final_s1_end[i+1];
				cluslen[len]=p1length1;
				switcarray[len]=switc+1;
				fgs+=1;
			}
		}

		//If segments do not merge
		else {
			switcarray[len]=switc+1;
			switc=0;

			for(k=1;k<all_kmer+1;k++) {
				hash1[hcounter]=freqa[k];
				hcounter+=1;
			}
	
			if(i+1==sg1) { //for last segment
				final_groups_start[fgs]=final_s1_start[i];
				final_groups_end[fgs]=final_s1_end[i];
				cluslen[len]=p1length1;

				//fprintf(output2,"final_groups_start[i] %d final_groups_end[i] %d\n",final_groups_start[fgs], final_groups_end[fgs] );
				len+=1; fgs+=1;
				final_groups_start[fgs]=final_s1_start[i+1];
				final_groups_end[fgs]=final_s1_end[i+1];
				cluslen[len]=final_s1_end[i+1]-final_s1_start[i+1]-2;
				switcarray[len]=1;
				//fprintf(output2,"fgsi %d final_groups_start[i] %d final_groups_end[i] %d\n",fgs,final_groups_start[fgs], final_groups_end[fgs] );
				total_cluster_length[fgs]=final_groups_end[fgs]-final_groups_start[fgs];
				fgs+=1;

				for(k=1; k<all_kmer+1; k++) {
					hash1[hcounter]=freqb[k];
					hcounter+=1;
				}
			}

			else { //for all other segments, save segments for non-contiguous clustering
				final_groups_start[fgs]=final_s1_start[i];
				final_groups_end[fgs]=final_s1_end[i];
				cluslen[len]=p1length1;
				total_cluster_length[fgs]=final_groups_end[fgs]-final_groups_start[fgs];
				len+=1; fgs+=1;
		    	}
			//fprintf(output2,"1final_groups_start[i] %d final_groups_end[i] %d\n",final_groups_start[fgs-1], final_groups_end[fgs-1] );
		}
	}

	sg1=1;
	for(i=1; i<fgs+1; i++) {
		final_s1_start [sg1]=final_groups_start[i];
		final_s1_end[sg1]=final_groups_end[i]; //fprintf(output2_2, "cluslen: %lf\tswitcarray: %d\n", cluslen[i],switcarray[i]);
		sg1 += 1; if (verb==1) { fprintf(output2, "final_groups_start[i] %d final_groups_end[i] %d\n",final_groups_start[i], final_groups_end[i]); }
	}
	count=0; count1=0; len=1; hcounter=1;
	//printf ("sg1 is %d\n",sg1);
	sg1=sg1-2; fgs=1;

	for (i=1; i<sg1; i++) {
        	if (switc==0){ p1length1= final_s1_end [i]-final_s1_start[i]+1-3; }

		p2length=final_s1_end [i+1]-final_s1_start[i+1]-1;
		np2length=p2length-1;
		total_length=p1length1+p2length;

 		for(k=1; k<all_kmer+1; k++) {
			if (switc==0) {
				for(s=final_s1_start [i]; s<final_s1_end[i]-2; s++) {
					if(k==hash[s]) { count=count+1; }
				}
				freqa[k]=count;
			}

			for(s=final_s1_start [i+1]; s<final_s1_end[i+1]-2; s++) {
				if(k==hash[s]) { count1=count1+1; }
			}
			freqb[k]=count1;
			freqab[k]=freqa[k]+freqb[k];
			count=0; count1=0;
		}
  
		enta = -((entropy(freqa,p1length1))/(log(2.0)));
		entb = -((entropy(freqb,p2length))/(log(2.0)));
		entab = -((entropy(freqab,total_length))/(log(2.0)));
		pi1 = (double)p1length1/total_length;
		pi2 = (double)1-pi1;
		jsd = entab-(pi1*enta)-(pi2*entb);

		if (m==2){ a=2.39; b=-7.66; c1=0.0029; d=0.841;}
		else if (m==0)	{ a=2.7784; b=-7.97084; c1=0.0; d=0.80; }				
		else if (m==1){ a=2.543; b=-4.77; c1=0.0; d=0.848; }
							
		beta=(c1*log(total_length)) + d;
		neff=(a*log(total_length)) + b;
		sarg=(log(2.0)*total_length*jsd);
		sx=gammp(dof/2,sarg);	
		smax=sx;
		//printf("%d %d %d %d %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %.11lf %lf %lf %lf %lf %.15lf\n",final_s1_start[i],final_s1_end [i],final_s1_start[i+1],final_s1_end[i+1],entab, enta, p1length1, pi1, pi1*enta, entb,p2length,pi2, pi2*entb,total_length, jsd,beta,neff,sarg,sx,smax);

		if (smax<thres1) {
			switc+=1;
			if(switc==1) {
				p1length1=(double)(total_length-1)/2;
				newp1l=total_length-1;

				for(k=1;k<all_kmer+1;k++) {
					freqa[k]=(double)((freqab[k])/2);
					tempafreq[k]=freqab[k];
				}
			}

			else {
				p1length1=(double)(newp1l+np2length)/(switc+1);
				newp1l=newp1l+p2length-1;

				for(k=1;k<all_kmer+1;k++) {
					freqa[k]=(double)(tempafreq[k]+freqb[k])/(switc+1);
					tempafreq[k]=tempafreq[k]+freqb[k];
				}
			}

			final_s1_start[i+1]=final_s1_start[i];

			if (i+1==sg1) {
				final_groups_start[fgs]=final_s1_start[i];
				final_groups_end[fgs]=final_s1_end[i+1];
				cluslen[len]=p1length1;
				switcarray[len]=switc+1;
				fgs+=1;
			}
		}

		else{
			switcarray[len]=switc+1;
			switc=0;

			for(k=1;k<all_kmer+1;k++){
				hash1[hcounter]=freqa[k];
				hcounter+=1;
	    		}

			if(i+1==sg1){	
				final_groups_start[fgs]=final_s1_start[i];
				final_groups_end[fgs]=final_s1_end[i];
				cluslen[len]=p1length1;

				//fprintf(output2,"final_groups_start[i] %d final_groups_end[i] %d\n",final_groups_start[fgs], final_groups_end[fgs] );
				len+=1; fgs+=1;
				final_groups_start[fgs]=final_s1_start[i+1];
				final_groups_end[fgs]=final_s1_end[i+1];
				cluslen[len]=final_s1_end[i+1]-final_s1_start[i+1]-2;
				switcarray[len]=1;
				//fprintf(output2,"fgsi %d final_groups_start[i] %d final_groups_end[i] %d\n",fgs,final_groups_start[fgs], final_groups_end[fgs] );
				total_cluster_length[fgs]=final_groups_end[fgs]-final_groups_start[fgs];
				fgs+=1;

				for(k=1;k<all_kmer+1;k++){
					hash1[hcounter]=freqb[k];
					hcounter+=1;
				}
			}

			else{
				final_groups_start[fgs]=final_s1_start[i];
				final_groups_end[fgs]=final_s1_end[i];
				cluslen[len]=p1length1;
				total_cluster_length[fgs]=final_groups_end[fgs]-final_groups_start[fgs];
				len+=1; fgs+=1;
			}
			if (verb==1) { fprintf(output2_1,"1final_groups_start[i] %d final_groups_end[i] %d\n",final_groups_start[fgs-1], final_groups_end[fgs-1] ); }
	    	}
	}

/********************************************************contiguous clustering ends*******************************************************/

	for(i=1; i<fgs; i++) {
		//fprintf(output2_3, "cluslen: %lf\tswitcarray: %d\n", cluslen[i],switcarray[i]);
		cluslen[i]=cluslen[i]*switcarray[i];
		//printf("cluslen: %lf\tswitcarray: %d\n", cluslen[i], switcarray[i]);
	}

/******************************************************non contiguous clustering starts***************************************************/
	count=0;
	int counter=1, temp2=0,v=1,jj=1,ii=1, nn,secondsegmentcounter=0,firstsegmentcounter=0,is=0,secondsegstart=1, cl=1, cluster [1000], seginfoclus [1000],cd=1,clusterdone[1000]={0},clusterdone1[1000]={0}, alonecluster=0,aaa=1,bbb=1,ccc=1, ascend1, switc1=0, switc2=0, tempswitcarrayis,tempcluslenis,switc2secondsegmentcounter,vv, tempq=9999999,tempclus1[1000],tmpcl1=1,cdremove,switc3,temp2q,tempis, clustera[5000], clusterb[5000], clacounter=1, newclusterb[5000], clusterid[5000], cid=1, cid1=1, nwclb=1, pos1=1, newclusterb1[5000], nwclb1=1,diff,clusterid1[5000];

	double p1length, p2length2,new1p1l, numerator, denominator;

	double *part1=malloc((mmk+1)*sizeof(double));
	double *part2=malloc((mmk+1)*sizeof(double));
	double *part12=malloc((mmk+1)*sizeof(double));
	double *temp1freq=malloc((mmk+1)*sizeof(double));
	double *nume=malloc(1000*sizeof(double));
	double *dene=malloc(1000*sizeof(double));
	double *temppart1freq=malloc((mmk+1)*sizeof(double));

	enta=0;entb=0;entab=0;p1length=0;p2length=0;np2length=0;total_length=0;firstsegmentcounter=0;

	for(is=1; is<fgs; is++) {//foreach segment-cluster

		for(i=1; i<fgs-1; i++) {// check if 1st segment is already clustered //is<fgs-1
			if(is==clusterdone[i]){
				v+=all_kmer;
				is+=1;
				firstsegmentcounter+=1;
			}
		}
	
		alonecluster=0;
		vv=(is*all_kmer)-all_kmer+1;
		v=(is*all_kmer)-all_kmer+1;

        	if(switc1==0) {
			for(nn=v; nn<v+all_kmer; nn++) { //for 1st segment get frequency
				part1[jj]=hash1[nn];
        			jj++;
			}

			p1length=cluslen[is]/switcarray[is]; //fprintf(output3, "%lf %lf %d\n",p1length,cluslen[is],switcarray[is]);
		}
			
		jj=1; secondsegstart=v;
		firstsegmentcounter+=1; secondsegmentcounter=is+1;

		for(q=secondsegmentcounter;q<fgs; q++) {//for second segment     
        		ii = 1; secondsegstart += 64;

			if(switc2==0) {

				for(i=1; i<fgs+1; i++) { //compare part2
                        				
					if(q==clusterdone1[i]) { //check if segment already clustered                       
                        			secondsegstart+=64;
						q+=1;
					}
				}
			}
 
			else{ q=tempq; }

			if(q<fgs) {
				ii=1;	
				secondsegstart=(q*all_kmer)-all_kmer+1;

    				for(j=secondsegstart; j<secondsegstart+64; j++) {//collect segemnt 2 frequency
        				part2[ii]=hash1[j];
        				ii+=1;		
				} 

        			for(i=1; i<all_kmer+1; i++) {//add segemnt 1 and 2 freq
           				part12[i]=part1[i]+part2[i];
				}

				p2length2=cluslen[q]/switcarray[q];	   								
				np2length=p2length2-1;
				total_length=p1length+p2length2;
				//fprintf(output3_1, "%lf  %lf  %lf  %lf  %d  %d  %lf\n",p1length,p2length2,cluslen[is],cluslen[q],switcarray[is],switcarray[q],total_length);

				enta = -((entropy(part1,p1length))/(log(2.0)));
				entb = -((entropy(part2,p2length2))/(log(2.0)));
				entab = -((entropy(part12,total_length))/(log(2.0)));
				pi1 = (double)p1length/total_length;
				pi2 = (double)1-pi1;
				jsd = entab-(pi1*enta)-(pi2*entb);
						
				if (m==2) { a=2.39; b=-7.66; c1=0.0029; d=0.841; }
				else if (m==0) { a=2.7784; b=-7.97084; c1=0.0; d=0.80; }				
				else if (m==1) { a=2.543; b=-4.77; c1=0.0; d=0.848; }

				beta=(c1*log(total_length)) + d;
                        	neff=(a*log(total_length)) + b;
		        	sarg=(log(2.0)*total_length*jsd);
                        	sx=gammp(dof/2,sarg);	
				smax=sx;
							
				//fprintf(output3_2, "%lf %lf %lf %d %d %.12lf %.18lf\n",total_length,p1length,p2length2,switcarray[is],switcarray[q], jsd,smax);

				//cluster if divergence is not significant
				if (smax<thres3) { 
      
					//printf("cluster %lf with %lf\n", p1length,p2length2);
					clustera[clacounter]=is; clusterb[clacounter]=q;
					clacounter+=1;

					if(q==tempq) { tempq=99999; }

					if(tmpcl1>1) {
						if(is==1) {
							for(i=1;i<tmpcl1;i++) {
								if(tempclus1[i]==q) {
									tempclus1[i]=0;
								}
							}
						}
					}
					switc1+=1;
					secondsegmentcounter=is+1;
					secondsegstart=v;
					if(switc1==1) {
						p1length=(double)((cluslen[is]+cluslen[q])/(switcarray[q]+switcarray[is]));		
						tempcluslenis=cluslen[is];						
						new1p1l=(cluslen[is])+cluslen[q];
						cluslen[is]=(cluslen[is])+cluslen[q];
							
						for(k=1; k<all_kmer+1; k++) {
							temppart1freq[k]=part1[k];
							numerator=(double)((part1[k]*(switcarray[is]))+(part2[k]*(switcarray[q])));
							denominator=(double)(switcarray[q]+switcarray[is]);
							temp1freq[k]=(double)((part1[k]*(switcarray[is]))+(part2[k]*(switcarray[q])));
							part1[k]=numerator/denominator;
							nume[k]=numerator/denominator;
						}
						tempswitcarrayis=switcarray[is];
						switcarray[is]=switcarray[is]+switcarray[q];
					}
					else {		
						p1length=(double)(cluslen[is]+cluslen[q])/(switcarray[q]+switcarray[is]);
						new1p1l=cluslen[is]+cluslen[q];
						cluslen[is]=cluslen[is]+cluslen[q];
							
						if(switc2==1) {	
							for(k=1; k<all_kmer+1; k++) {
								temp1freq[k]=(double)(part1[k]*switcarray[is]);
							}
						}
								
						for(k=1; k<all_kmer+1; k++) {
							temppart1freq[k]=part1[k];
							part1[k]=(double)(temp1freq[k]+(part2[k]*(switcarray[q])))/(switcarray[q]+switcarray[is]);
							temp1freq[k]=(temp1freq[k]+(part2[k]*(switcarray[q])));
						}
						tempswitcarrayis=switcarray[is];
						switcarray[is]=switcarray[is]+switcarray[q];
						switc2=0;
					}//else ends

					alonecluster += 1;		
					clusterdone[q]=q;clusterdone1[q]=q;
					cd += 1; bbb += 1; aaa += 1;
					secondsegmentcounter=is+1;

					if(firstsegmentcounter>1) {		
						switc2=1;
						switc2secondsegmentcounter=is;
						j=(is*all_kmer)-all_kmer+1;
						clusterdone[q]=q; 
						clusterdone1[q]=q;	
						tempclus1[tmpcl1]=is;
						tmpcl1+=1;

						if(is<tempq) { tempq=is; }

						for(k=1; k<all_kmer+1; k++) {
							temp1freq[k]=(double)((part1[k]*(switcarray[is]))+(part2[k]*(switcarray[q])));
							hash1[j]=part1[k];
							j++;
						}

						j=(is*all_kmer)-all_kmer+1;
						p2length2=(double)(cluslen[is]+cluslen[q])/(switcarray[is]+switcarray[q]);

						for(k=1; k<all_kmer+1; k++) {		
							part2[k]=part1[k];
						}

						tempis=is;
						firstsegmentcounter=1; is=1;
						p1length=(double)(cluslen[is]/switcarray[is]);
								
						for(k=1; k<all_kmer+1; k++){
							part1[k]=hash1[k];//part1[k]=(dene[k]);
						}
					}// if firstsegmentcounter close
					q=1;		
	    			}//if cluster close

				//If segments do not merge in same cluster
				else {
									
					switc1=0;
					clusterdone[is]=0;
								
					if(switc2==1) {
								
						secondsegstart=1;
						clusterdone[q]=0;	
						j=(switc2secondsegmentcounter*all_kmer)-all_kmer+1;
								
						if(tempis==1) {	
							for(i=1;i<tmpcl1;i++) {
								if(tempclus1[i]>q) {
									cdremove=tempclus1[i];
									clusterdone1[cdremove]=0;
								}
							}
						}

						clusterdone1[q]=0; switc3=1;
						if(switc3==1) { q=q; //printf("switc3 %d q %d\n", switc3,q); 
						}
						else { q=1; switc3=0; }						
					}		
					switc2=0;			
				}
			 }
		}//for second segment close

		j=(is*all_kmer)-all_kmer+1;

		for(k=1; k<all_kmer+1; k++) {
			dene[k]=nume[k];
			hash1[j]=part1[k];
			j++;
		}

		cluslen[is]=(double)(p1length*switcarray[is]);
		aaa+=1;
		bbb+=1;
		v+=64;
		cl++;
	}//foreach segment close

/**************************************************non-contiguous clustering ends*********************************************************/

	//Assign cluster ids to segments based on mergers
	for(i=1; i<clacounter; i++) {
		for(j=1; j<clacounter; j++) {
			if(clustera[i]==clusterb[j]) {
				clustera[i]=clustera[j];
			}
		}
	}

	for(i=1; i<clacounter; i++) {
		for(j=i+1; j<clacounter; j++) {
			if(clustera[i]>clustera[j]) {
				ascend= clustera[i];
				ascend1=clusterb[i];
				clusterb[i]=clusterb[j];
				clustera[i]=clustera[j];
				clustera[j]=ascend;
				clusterb[j]=ascend1;
			}
		}
		//printf("ascend clustera [i] %d cluster b[i] %d\n", clustera[i], clusterb[i]);
	}

	pos1=1;

	for(i=1; i<clacounter; i++) {
		if(clustera[i]<clustera[i+1]) {		
			for(j=pos1; j<i+1; j++) {
				newclusterb[nwclb]=clusterb[j];
				clusterid[cid]=cid1;
                                cid += 1; nwclb += 1;
			}

			newclusterb[nwclb]=clustera[i];
			clusterid[cid]=cid1;
			cid1+=1; cid+=1;
			nwclb+=1; pos1=i+1;		
		}
	}

	for(i=pos1; i<clacounter; i++) {
		newclusterb[nwclb]=clusterb[i];
        	nwclb+=1;
		clusterid[cid]=cid1;
		cid+=1;
	}
		
	newclusterb[nwclb]=clustera[clacounter-1];
	clusterid[cid]=cid1;
	cid+=1;cid1+=1;

	//sort ascending
	for(i=1;i<cid-1;i++){
		for (j = i+1; j<cid; j++) {
			if (newclusterb[i]>newclusterb[j]) {
				ascend= newclusterb[i];
				ascend1=clusterid[i];
				newclusterb[i] = newclusterb[j];
				clusterid[i] = clusterid[j];
               			newclusterb[j] = ascend;
				clusterid[j] = ascend1;
			}
		}
		//printf(" final output %d %d \n", clusterid[i],newclusterb[i] );
	}

	pos1=1;

	for (i=1; i<fgs;i++){
		for (j=1; j<cid;j++){
			if (i==newclusterb[j]){
				clusterid1[i]=clusterid[j];
				//printf("newclusterb[j] %d clusterid1[%d] %d clusterid[%d] %d\n",newclusterb[j],i,clusterid1[i],j, clusterid[j]);
			}
		}
	}

	for (i=1; i<fgs;i++){
		if (clusterid1[i]==0){
			clusterid1[i]=cid1;
			cid1++;
		}
	}

	//Print segments with cluster ids
	FILE *output4 = fopen("segmentation_clustering.txt", "w"); 
	for(i=1;i<fgs;i++){
		fprintf(output4,"%d\t%d\t%d\n",final_groups_start[i], final_groups_end[i],clusterid1[i] ); 
	} 

	fclose(output1);
	fclose(output2_1);
	free(part1);
	free(part2);
	free(part12);
	free(tempafreq);
	free(temp1freq);
	free(cluslen);
	free(freqab);
	free(freqa);
	free(freqb);
	free(total_cluster_length);
	free(nume);
	free(dene);
	free(temppart1freq);
	//fclose(output3);
	fclose(output4);
}


/************************************************************segmentation starts**********************************************************/

int segmentation(int *hash, int parent_start, int parent_end) {
	int i=0, j=0, k, q, length,  seg_coord=0, s1_start, length1, p2length, p1length=0, s1_end, inter, newparentstart, newparentend, greater; 
	double ent, ent1, ent2, pi1, pi2,  n=0, maxi=-999999999.0, invlen1;
	double a, b, c, d, beta, neff, arg, sx, smax, sarg;
	double distance1;

	double *freq = malloc((mmk+1)*sizeof(double));
	double *freq1 = malloc((mmk+1)*sizeof(double));
	double *freq2 = malloc((mmk+1)*sizeof(double));
	
	for (k=1; k<all_kmer+1; k++) {
		freq1[k] = 0; freq2[k] = 0;
		}

	length = parent_end-parent_start+1;
	inter = length/10000;
	if(inter < 1) { inter = 1; }
	newparentstart = parent_start+mmk; 
	newparentend = parent_end-mmk;
	p1length = newparentstart-parent_start+1;
	p2length = length-p1length-4;
	length1 = p1length+p2length;
	invlen1 = (double)1/(length1);
	
	for (i=newparentstart; i<newparentend; i+=inter) {
		for (k=1; k<all_kmer+1; k++) {
			freq1[k] = 0; freq2[k] = 0;
		}
		//Count frequencies of trinucloetides
		for (j=parent_start; j<i-2; j++) {
			freq1[hash[j]] += 1;
		}

		for (q=i+1; q<parent_end-2; q++) {
			freq2[hash[q]] += 1;
		}

		for (k=1; k<all_kmer+1; k++) {
			freq[k] = freq1[k]+freq2[k];
		}

		//Calculate entropy
		ent = -((entropy(freq, length1))/(log(2.0)));
		ent1 = -((entropy(freq1, p1length))/(log(2.0)));
		ent2 = -((entropy(freq2, p2length))/(log(2.0)));
		pi1 = (double)p1length*invlen1;
		pi2 = 1-pi1;
		//Calculate divergence between the segments
		distance1 = ent-(pi1*ent1)-(pi2*ent2);

		//Get position of maximum divergence
		if(maxi < distance1) {
			maxi = distance1;
			seg_coord = i;
		}
		
		p1length += inter;
		p2length -= inter;
	}
	
	s1_start = parent_start;
	s1_end = seg_coord;
	
	//Parameters for Markov models
	if (m == 2) { a=2.39; b=-7.66; c=0.0029; d=0.841; }
	else if (m == 0) { a=2.7784; b=-7.97084; c=0.0; d=0.80; }	         						
	else if (m == 1) { a=2.543; b=-4.77; c=0.0; d=0.848; }

	//Do significance testing
	beta = (c*log(length)) + d;
	neff = (a*log(length)) + b;
	sarg = (log(2.0)*beta*length*maxi);
	sx = gammp(dof/2,sarg);	
	smax = pow(sx,neff);

	//---------------------------------------------------------------if significant----------------------------------------------------
	//Continue segmentation
	if(smax > thres) {
		s2_start[s2j] = s1_end+1;
		s2_end[s2j] = parent_end;
		s2j += 1;
		segmentation(hash, s1_start, s1_end);
	}

	//-------------------------------------------------------------if not significant--------------------------------------------------
	//Do not segment further
	else {
		final_s1_start[sg1] = parent_start;
		final_s1_end[sg1] = parent_end; 
		//printf("from function final_sement_one%d start:%d end:%d\n",sg1, final_s1_start[sg1], final_s1_end[sg1]);
		if (verb == 1) { printf("%d\t%d\t%.12lf\n", parent_start, parent_end, distance1); }
		sg1 += 1;
		//Get next segment for segmentation
		for(k1=l; k1<s2j; k1++) {
	   		l += 1;
	   		segmentation(hash, s2_start[k1], s2_end[k1]);
		}	
	}

	free(freq); free(freq1); free(freq2); 
}

/************************************************************* Entropy *******************************************************************/
//Function to calculate entropy. Takes input as frequency array

double entropy(double *freq_kmer, double length) {
	int i=1,i1=1;
	double ent=0.0,sumkmer=0.0, ent1=0.0,wordprob, *probkmer=malloc((mmk+1)*sizeof(double));

	for (i=1;i<all_kmer+1;i++) {
                sumkmer=(freq_kmer[i1]+freq_kmer[i1+1]+freq_kmer[i1+2]+freq_kmer[i1+3]); //Get sum of dinucleotides (kmer-1)	
 		if (sumkmer>0) {
			probkmer[i]=(double)freq_kmer[i]/sumkmer; //Calculate probability of trinucleotides (kmers)
                     	if (probkmer[i]>0.0) {
				ent1=ent1+probkmer[i]*((log (probkmer[i])));
			}	
 		}
                
		if (i%4 ==0) {
			wordprob=(double)sumkmer/length; //calculate probability of dinucleotide (kmer-1)
			ent+=ent1*(wordprob);  
			ent1=0.0;
			i1=i1+4;
		}
	}
	free(probkmer);
	return ent;
}
/*****************************************************************************************************************************************/
float gammp(float a, float x) {
	void gcf(float *gammcf, float a, float x, float *gln);
	void gser(float *gamser, float a, float x, float *gln);
	void nrerror(char error_text[]);
	float gamser, gammcf,gln;

//	if(x < 0.0 || a <= 0.0) nrerror("invalid arguments in routine gammp ");
//	if(x < 0.0 || a <= 0.0) nrerror("");
	if (x< (a+1.0)) {
		gser(&gamser, a,x,&gln);
		return gamser;
		}
	else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}

/*****************************************************************************************************************************************/
void gser(float *gamser, float a,float x, float *gln) {
	//float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;
	*gln=gammln(a);
	if (x <= 0.0) {
		*gamser=0.0;
		return;
	}
	else {
		ap=a;
		del=sum=1.0/a;
		for(n=1;n<=ITMAX;n++) {
			++ap;
			del*=x/ap;
			sum+=del;
			if(fabs(del)<fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		return;
	}
}

/*******************************************************************************************************************************************************************************************************/
void gcf(float *gammcf, float a, float x, float *gln) {
	//float gammln (float xx);
	void nrerror (char error_text[]);
	int i;
	float an,b,c,d,del,h;
	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for(i=1;i<=ITMAX;i++) {
		an=-i*(i-a);
		b+=2.0;
		d=an*d+b;
		if(fabs(d)<FPMIN) { d=FPMIN; }
		c=b+an/c;
		if(fabs(c)< FPMIN) { c=FPMIN; }
		d=1.0/d;
		del=d*c;
		h*=del;
		if (fabs(del-1.0)<EPS) { break; }
	}

	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}

/*****************************************************************************************************************************************/
float gammln(float xx) {
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091, -1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp-=(x+0.5)*log(tmp);
	ser=1.000000000190015;
	for(j=0;j<=5;j++) {
		ser+= cof[j]/++y;
	}
	return -tmp+log(2.5066282746310005*ser/x);
}
/*****************************************************************************************************************************************/
void nrerror (char error_text[])
	{
	printf("%s\n",error_text);
	}
/*****************************************************************************************************************************************/

