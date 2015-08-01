#!/usr/bin/env python
# -*- coding: ASCII -*-

#####This python script that takes a list of bamfiles and performs variant calling such that the output is in the format required for lcMLkin to calculate pairwise relatedness.
#####Please note it is very slow, and is only designed for calling genotypes for specific SNPs. I do not recommend using this as a general variant caller.
#####We are working on a better parser (this was an experimental one developed to test lcMLkin), but for now, if using 100 samples with 100K loci, it may need to run overnight.
#####Pysam and numpy must be installed in the version of python used
#####Bam files must be indexed. The SNPfile must have the tab seperated fields in the following order: chromosome, position (one-based), reference allele, alternate allele
#####Written (poorly) by Krishna Veeramah

#####usage is ./SNPbam2vcf.py <bamlist> <fileoutname> <target_SNP_file>'


###import libraries
import string
import numpy as np
import pysam
import gzip
import math
import copy
from sys import argv


###Input arguments
BAMin=argv[1] #flat text file consisiting of list of bams (must be indexed), one per line
filenameout=argv[2] #creates a vcf
SNPfile=argv[3] #must have the tab seperated fields in the following order chromosome, position (one-based), reference allele, alternate allele
BC=0 #choose whether to use the bayesian caller for individual genotypes (0=no, 1=yes). I don't recommend turning this on the current  version, as it's quite slow, and allele frequencies shouldn't really be biased if the genotypes are incorreclty called, especially as the genotype calling is not reference aware
MQ_t=20 #mapping quality threshold
BQ_t=5 #base_qualitythreshold
GQ_t=0.1 #GQ threshold
min_sams=1 #minimum samples with reads at locus needed to attempt bayesian calling of genotypes


###converts phred score to probability
def phred2prob(x):
    return 10.0**(-x/10.0)

###converts probability to phred score
def prob2phred(x):
    return -10*math.log10(x)

###diploid caller assuming all alternative alleles are possible (error is divided by three for now, could add a more complex model)
def geno_caller_10GT(X):
   
    GL=np.zeros(10)   #all 10 possible genotypes and order = AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
    hap=np.zeros((len(X),4))  #all 4 haploid possibilities, A,C,G,T

    all_dic={}
    all_dic['A']=0
    all_dic['C']=1
    all_dic['G']=2
    all_dic['T']=3
    
    count=0
    for g in range(len(X)):
        if all_dic.has_key(X[g][0])==False:
            continue
        hap[g]=phred2prob(X[g][1])*(1.0/3.0)
        hap[g][all_dic[X[g][0]]]=1-phred2prob(X[g][1])

        GL[0]=GL[0]+math.log10(hap[g][0])
        GL[1]=GL[1]+math.log10((hap[g][0]+hap[g][1])/2)
        GL[2]=GL[2]+math.log10((hap[g][0]+hap[g][2])/2)
        GL[3]=GL[3]+math.log10((hap[g][0]+hap[g][3])/2)

        GL[4]=GL[4]+math.log10(hap[g][1])
        GL[5]=GL[5]+math.log10((hap[g][1]+hap[g][2])/2)
        GL[6]=GL[6]+math.log10((hap[g][1]+hap[g][3])/2)

        GL[7]=GL[7]+math.log10(hap[g][2])
        GL[8]=GL[8]+math.log10((hap[g][2]+hap[g][3])/2)

        GL[9]=GL[9]+math.log10(hap[g][3])
        count+=1

    if count==0:
        GL.fill(-9)
    return GL


###extract reads for a give position in a bam
def extract_bam_SNP(samfile,chromo,pos,BQ,MQ):
    var_list=[]
    for pileupcolumn in samfile.pileup(chromo,pos-1,pos):
        if pileupcolumn.pos==pos-1:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    if pileupread.alignment.mapping_quality>=MQ:
                        if pileupread.alignment.is_duplicate==False:
                            if ord(pileupread.alignment.qual[pileupread.query_position])-33>=BQ:
                                var_list.append([pileupread.alignment.query_sequence[pileupread.query_position],ord(pileupread.alignment.qual[pileupread.query_position])-33])
    return var_list


###Version of Depristo's bayesian caller
def bayes_caller(X):
    het=0.001
    nb_ind=len(X)


    thetas=np.zeros(((nb_ind*2)+1),dtype='float32')  #matrix with theta priors
    for g in range(1,len(thetas)):
        thetas[g]=het/float(g)
    thetas[0]=1-sum(thetas[1:])

    best_genos=np.zeros(((nb_ind*2)+1,nb_ind),dtype='int64')  #matrix for storing best genotypes for a given allele frquency
    best_genos_ll=np.zeros(((nb_ind*2)+1,nb_ind),dtype='float32')  #matrix for storing best genotypes for a given allele frquency

    p_best_genos=np.zeros(((nb_ind*2)+1),dtype='float32')  #matrix for storing the likelihood for each p

    for g in range(len(X)):   ### p = 0
        best_genos_ll[0][g]=best_genos_ll[0][g]+X[g][0]

    p_best_genos[0]=np.sum(best_genos_ll[0])

    ##perform greedy algorithm to find best likelihoods for each possible configuration as described in Depristo
    for g in range(1,len(p_best_genos)):
        p_best_genos[g]=-1000000000000
        genos=np.copy(best_genos[g-1])
        genos_ll=np.copy(best_genos_ll[g-1])
        for gg in range(len(genos)):
            if genos[gg]<2:
                genos_mod=np.copy(genos)
                genos_mod[gg]+=1
                genos_ll_mod=np.copy(genos_ll)
                genos_ll_mod[gg]=X[gg][genos[gg]+1]
                ll=np.sum(genos_ll_mod)
                if ll>p_best_genos[g]:
                    p_best_genos[g]=ll
                    best_genos[g]=np.copy(genos_mod)
                    best_genos_ll[g]=np.copy(genos_ll_mod)


    ##Calculate variant quality
    #denom=np.max(p_best_genos)+np.log10(np.sum(10**(p_best_genos-np.max(p_best_genos))))
    #phred=-10*(np.log10(thetas[0])+p_best_genos[0]-denom)
    marg_comp=np.log10(thetas)+p_best_genos
    denom=np.max(marg_comp)+np.log10(np.sum(10**(marg_comp-np.max(marg_comp))))
    phred=-10*(np.log10(thetas[0])+p_best_genos[0]-denom)

    VQ=int(round(phred))
    
    geno_comb=best_genos[np.argmax(p_best_genos)]
    
    if np.sum(geno_comb)==0:
        stat='REF'

    else:
        stat='ALT'
    return geno_comb,VQ,stat

###dictionary for GL order
geno_ord={}
geno_ord['AA']=0
geno_ord['AC']=1
geno_ord['AG']=2
geno_ord['AT']=3
geno_ord['CA']=1
geno_ord['CC']=4
geno_ord['CG']=5
geno_ord['CT']=6
geno_ord['GA']=2
geno_ord['GC']=5
geno_ord['GG']=7
geno_ord['GT']=8
geno_ord['TA']=3
geno_ord['TC']=6
geno_ord['TG']=8
geno_ord['TT']=9

###dictionary for vcf genotype codes
vcf_gt={}
vcf_gt[0]='0/0'
vcf_gt[1]='0/1'
vcf_gt[2]='1/1'


###Read bamlist
file=open(BAMin,'r')
samples=file.read()
samples=string.split(samples,'\n')
if samples[-1]=='':
    del(samples[-1])


###Read SNPlist
file=open(SNPfile,'r')
SNPs=file.read()
SNPs=string.split(SNPs,'\n')
if SNPs[-1]=='':
    del(SNPs[-1])
   

sample_name=[]
###Get sample name for each bam. Each bam must only contain one sample, otherwise the program will take the sample in first readgroup in the header
for g in range(len(samples)):
    samfile = pysam.AlignmentFile(samples[g], "rb")
    sample_name.append(samfile.header['RG'][0]['ID'])

###start writing outfile with header
outfile=open(filenameout,'w')
out='##Custom variant caller for_lcMLkin\n##Bayesian_caller='+str(BC)+',MappingQuality_filter='+str(MQ_t)+',BaseQuality_filter='+str(BQ_t)+',GenotypeQuality_filter='+str(GQ_t)+',bamlist='+BAMin+',filenameout='+filenameout+',SNPfile='+SNPfile+'\n'
out=out+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
for x in sample_name:
    out=out+'\t'+x
    
out=out+'\n'

outfile.write(out)



###Iterate through SNPs
for g in range(len(SNPs)):
    k=string.split(SNPs[g])
    chromo=k[0]
    pos=int(k[1])
    ref=k[2]
    alt=k[3]

    ###iterate through samples
    GLs=[]
    GLs_nomiss=[]
    GQs=[]
    GTs=[]
    GTs=[]
    AD=[]
    for gg in range(len(samples)):
        ###extract reads for a position
        samfile = pysam.AlignmentFile(samples[gg], "rb")
        read_list=extract_bam_SNP(samfile,chromo,pos,BQ_t,MQ_t)

        ###get allele depth infomation
        if len(read_list)>0:
            alls=list(zip(*read_list))[0]
            AD.append([alls.count(ref),alls.count(alt)])
        else:
            AD.append([0,0])

        ###get the genotype likelihood for the 3 possible genotypes based on the ref and alt alleles provided
        ll=geno_caller_10GT(read_list)        
        GL=[ll[geno_ord[ref+ref]],ll[geno_ord[ref+alt]],ll[geno_ord[alt+alt]]]
        GLs.append(GL)
        
        ###calcule genotype quality
        GL_sort=copy.deepcopy(GL)
        GL_sort.sort()
        GQ=-10*math.log10(10**GL_sort[1]/10**GL_sort[2])
        GQs.append(GQ)

        ###assign genotype using highest likelihood
        if BC==1:
            if GQ>GQ_t:
                GLs_nomiss.append(GL)
            
        elif BC==0:
            if GQ>GQ_t:       
                GTs.append(vcf_gt[np.argmax(np.asarray(GL))])
            else:
                GTs.append('./.')

    if BC==1:
        if len(GLs_nomiss)>min_sams:  #must have genotype data for at least 20 individuals to consider SNP
            genos_b,VQ,s=bayes_caller(GLs_nomiss)
            AF=(np.sum(genos_b)/(float(len(genos_b))*2))
            it=0
            for ggg in range(len(GLs)):
                if -9.0 in GLs[gg]:
                    GTs.append('./.')
                else:
                    GTs.append(vcf_gt[genos_b[it]])
                    it+=1
                
        else:
            AF=0.0
                
    ###estimate allele frequencies for single-base calling
    elif BC==0:
         AF=(GTs.count('0/1')+GTs.count('1/1')*2) / (float((len(GTs)-GTs.count('./.')))*2)
         VQ=100
         s='ALT'

    ###only write to file if site is variable (otherwise useless for kinship estimation)
    if 0.0<AF<1.0:                        
        out=k[0]+'\t'+k[1]+'\t.\t'+k[2]+'\t'+k[3]+'\t'+str(VQ)+'\t'+s+'\tALT_AF='+str(AF)+'\tGT:DP:AD:GL:GQ'

        for ggg in range(len(GTs)):
            out=out+'\t'+GTs[ggg]+':'+str(sum(AD[ggg]))+':'+str(AD[ggg][0])+','+str(AD[ggg][1])+':'+str(10**GLs[ggg][0])+','+str(10**GLs[ggg][1])+','+str(10**GLs[ggg][2])+':'+str(int(round(GQs[ggg])))

        out=out+'\n'
        outfile.write(out)                                                                                  

    if g%100==0:
        print 'Processed up to ' +chromo+':'+str(pos)+', SNP nb = '+str(g)

outfile.close()
