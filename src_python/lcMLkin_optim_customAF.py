##import time
import subprocess
from subprocess import Popen
import os
##import dircache
import string
from string import join
import math
##from time import strftime
import random
from random import randint
from random import uniform
from random import gauss
from random import gammavariate
from random import betavariate
from math import sqrt
from sys import argv
import numpy as np
import datetime
import numpy.ma as ma
from scipy.optimize import minimize
from scipy.optimize import fmin
import gzip

def kin(k,IBS,mask_matrix):
    k3=np.array([1-(k[0]+k[1]),k[0],k[1]])
    ll=-np.sum(ma.masked_array(np.log(np.sum(IBS*k3,axis=1)),mask=mask_matrix))
    pen=0
    if k3[0]<0:
        pen+=1
    if k3[0]>1:
        pen+=1
    if k3[1]<0:
        pen+=1
    if k3[1]>1:
        pen+=1
    if k3[2]<0:
        pen+=1
    if k3[2]>1:
        pen+=1
    if 4*k3[2]*k3[0]>k3[1]**2:
        pen+=1
    if np.isinf(ll)==True:
        pen+=1
    if pen>0:
        ll=10E10
    return ll

def GLkin(k,GL,IBS,mask_matrix):
    k3=np.array([1-(k[0]+k[1]),k[0],k[1]])
    ll=-np.sum(ma.masked_array(np.log(GL[0]*np.sum(IBS[0]*k3,axis=1)+GL[1]*np.sum(IBS[1]*k3,axis=1)+GL[2]*np.sum(IBS[2]*k3,axis=1)+GL[3]*np.sum(IBS[3]*k3,axis=1)+GL[4]*np.sum(IBS[4]*k3,axis=1)+GL[5]*np.sum(IBS[5]*k3,axis=1)+GL[6]*np.sum(IBS[6]*k3,axis=1)+GL[7]*np.sum(IBS[7]*k3,axis=1)+GL[8]*np.sum(IBS[8]*k3,axis=1)),mask=mask_mat))
    pen=0
    if k3[0]<0:
        pen+=1
    if k3[0]>1:
        pen+=1
    if k3[1]<0:
        pen+=1
    if k3[1]>1:
        pen+=1
    if k3[2]<0:
        pen+=1
    if k3[2]>1:
        pen+=1
    if 4*k3[2]*k3[0]>k3[1]**2:
        pen+=1
    if np.isinf(ll)==True:
        pen+=1
    if pen>0:
        ll=10E10
    return ll

filenamein=argv[1]#'sim80_10K.vcf'
if filenamein[-2:]=='gz':
    file = gzip.open(filenamein)    
else:
    file = open(filenamein)
data=file.read()
data=string.split(data,'\n')
file.close()

if data[-1]=='':
    del(data[-1])

##unrels=['1_0', '1_1', '1_2', '1_3', '1_4', '1_5', '1_6', '1_7',
##        '2_0', '2_1', '2_2', '2_3', '2_4', '2_5', '2_6', '2_7',
##        '3_0', '3_1', '3_2', '3_3', '3_4', '3_5', '3_6', '3_7',
##        '4_0', '4_1', '4_2', '4_3', '4_4', '4_5', '4_6', '4_7',
##        '5_0', '5_1', '5_2', '5_3', '5_4', '5_5', '5_6', '5_7',
##        '6_0', '6_1', '6_2', '6_3', '6_4', '6_5', '6_6', '6_7',
##        '7_0', '7_1', '7_2', '7_3', '7_4', '7_5', '7_6', '7_7',
##        '8_0', '8_1', '8_2', '8_3', '8_4', '8_5', '8_6', '8_7',
##        '9_0', '9_1', '9_2', '9_3', '9_4', '9_5', '9_6', '9_7',
##        '10_0', '10_1', '10_2', '10_3', '10_4', '10_5', '10_6', '10_7']


head=string.split(data[0],'\t')[9:]
unrel_ind=[]

for g in range(len(head)):
#    if head[g] in unrels:
    unrel_ind.append(g)


for g in range(len(data)):
    if data[g][0]<>'#':
        start=g
        break


#########Find out order of INFO field
datum=[]
alls=[]####
for g in range(start,len(data)):
    k=string.split(data[g],'\t')
    datum.append(k[9:])
    alls.append(k[:9])####

GT_key=-9
GQ_key=-9
GL_key=-9

INFO=string.split(data[start],'\t')[8]
INFO=string.split(INFO,':')

for g in range(len(INFO)):
    if INFO[g]=='GT':
        GT_key=g
    if INFO[g]=='GQ':
        GQ_key=g
    if INFO[g]=='GL':
        GL_key=g

#########End of find out order of INFO field

data=[] #emtpy memory of initial file

nbSNPs=len(datum)

AF=np.zeros(nbSNPs,float)

##print 'calculating allele frequencies'
##for g in range(len(datum)):
##    counts=[]
##    for gg in range(len(unrel_ind)):
##        geno=string.split(datum[g][unrel_ind[gg]],':')[GT_key]
##        if geno <> './.':
##            counts.append(geno[0])
##            counts.append(geno[-1])
##    AF[g]=counts.count('0')/float(len(counts))

#####New code for giving a file with precomputed AFs'
filenameinSNP=argv[2]#'SNP_CDX.freq'
fileSNP = open(filenameinSNP)
data=fileSNP.read()
data=string.split(data,'\n')
fileSNP.close()

if data[-1]=='':
    del(data[-1])

if len(data)<>len(datum):
    print 'Size of SNP list doesnt match number of SNPs in vcf'
    huh
    
for g in range(len(data)):
    k=string.split(data[g],'\t')
    AF[g]=1-float(k[4])

    
print 'calculating prestored IBS|IBD'

##For every SNP (as each has a different allele frequency) calculate the P(IBS=x|IBD=z) for all combinations. We only do it for the 3 IBD values, i.e. assuming no inbreeding in the two individuals

##Makes a matrix to store these values. One matrix for each possible genotype combination we might observe. If we consider the reference as always p, there are 9 possible combinations.
##In a standard table like Mulligan or Thompson, some of these are collapsed (i.e. ppqq is the same as qqpp, as p has no meaning  has no meaning in that context)

#IBS0=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
ppqq=np.zeros((nbSNPs,3),float)
qqpp=np.zeros((nbSNPs,3),float)

#IBS1=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
pppq=np.zeros((nbSNPs,3),float)
pqpp=np.zeros((nbSNPs,3),float)
pqqq=np.zeros((nbSNPs,3),float)
qqpq=np.zeros((nbSNPs,3),float)

#IBS2=np.zeros((nbSNPs,3),float)   #IBD0,IBD1,IBD2
pppp=np.zeros((nbSNPs,3),float)
pqpq=np.zeros((nbSNPs,3),float)
qqqq=np.zeros((nbSNPs,3),float)


#populates matrix with actual probabilities
for g in range(len(AF)):
    p=AF[g]
    q=1.0-p
    #IBS0[g]=[2.0*(p**2)*(q**2),0,0]
    ppqq[g]=[(p**2)*(q**2),0,0]   
    qqpp[g]=[(q**2)*(p**2),0,0]
    
    #IBS1[g]=[4*(p**3)*q+4*p*(q**3),2*(p**2)*q+2*p*(q**2),0]
    pppq[g]=[(p**2)*(2*p*q),(p**2)*q,0]
    pqpp[g]=[(2*p*q)*(p**2),(p*q)*p,0]
    pqqq[g]=[(2*p*q)*(q**2),(p*q)*q,0]
    qqpq[g]=[(q**2)*(2*p*q),(q**2)*p,0]

    
    #IBS2[g]=[p**4+q**4+4*(p**2)*(q**2),p**3+q**3+(p**2)*q+p*(q**2),1]
    pppp[g]=[(p**4),p**3,p**2]
    pqpq[g]=[(2*p*q)*(2*p*q),p*q*(p+q),2*(p*q)]
    qqqq[g]=[(q**4),q**3,q**2]


#Store in a single matrix, one dimension for each of the 9 possible genotype combinations we might observe
    
IBS_all=np.array([ppqq,qqpp,pppq,pqpp,pqqq,qqpq,pppp,pqpq,qqqq])  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT



#pairwise analysis

k_combs=[]   #define parameter space i.e. what k1 and k2 coefficients (k0 is whatever is left over) we will look at. Currently set up as a grid moving 1% at a time, though bound by certain impossible constraints

k1_lis=[]
for g in range(101):
    k1_lis.append(g/100.0)

k2_lis=[]
for g in range(101):
    k2_lis.append(g/100.0)


for g in range(len(k1_lis)):
    k1=k1_lis[g]
    for gg in range(len(k2_lis)):
        k2=k2_lis[gg]       
        k0=(1.0-(k1+k2))
        if k1+k2<=1:
            if k0+k1+k2==1.0:
                if 4*k2*k0<k1**2:
                    k_combs.append([k0,k1,k2])

k_combs.sort()
k_combs=np.array(k_combs)

pw=[]   ##here we do every combinationpw
for g in range(len(head)):
    for gg in range(g+1,len(head)):
#        if string.split(head[g],'_')[0]==string.split(head[gg],'_')[0]:
        pw.append([g,gg])
                   

filenameout=filenamein+'.relateF_optim'
file=open(filenameout,'w') 
out='Ind1\tInd2\tZ0bg\tZ1bg\tZ2bg\tPI_HATbg\tZ0ag\tZ1ag\tZ2ag\tPI_HATag\tnbSNP\n'
file.write(out)
file.close()



print 'starting pairwise IBD computations\n'

print out[:-1]

#iterate through each pairwise comparison
for g in range(len(pw)):
    PIBS=np.zeros((nbSNPs,9),float)   ###matrix for all possible pairs of genotypes for every SNP. This will eventually store genotype likelihoods based on the calling likelihood (for example based on read depth)
    PIBSg=np.zeros((nbSNPs,9),float)   ###matrix all possible pairs of genotypes for every SNP. This will eventually store genotype likelihoods based on the calling likelihood (for example based on read depth

    ###A matrix to denote SNPs we may want to mask for two reasons(see below)
    mask_mat=np.zeros(nbSNPs,int)

    ##iterate through each SNP
    for gg in range(nbSNPs):
        ###mask SNPs where allele frequency is fixed
        if AF[gg]==1.0:
            mask_mat[gg]=1
        if AF[gg]==0.0:
            mask_mat[gg]=1
            
        l1=string.split(string.split(datum[gg][pw[g][0]],':')[GL_key],',')  ###pulls out info field for first individual from VCF, pulls out the three precomputed genotype likelihoods
        l2=string.split(string.split(datum[gg][pw[g][1]],':')[GL_key],',')  ###pulls out info field for second individual from VCF, pulls out the three precomputed genotype likelihoods

        ###converts likelihoods from strings to floating point numbers
        for ggg in range(len(l1)):
            l1[ggg]=float(l1[ggg])
            l2[ggg]=float(l2[ggg])
            if l1[ggg]==-9.0:   #if one of those likelihoods comes out as negative (some anomaly) we mask that SNPs
                mask_mat[gg]=1
            if l2[ggg]==-9.0:
                mask_mat[gg]=1

        ###We now calculate the probability of observing all possible two genotype combination by multiplying their likelihoods together
        #I0=(l1[0]*l2[2])+(l1[2]*l2[0])  #AA,TT    TT,AA
        PPQQ=(l1[0]*l2[2])
        QQPP=(l1[2]*l2[0])

        #I1=(l1[0]*l2[1])+(l1[1]*l2[0])+(l1[1]*l2[2])+(l1[2]*l2[1])       #AA,AT    AT,AA   AT,TT   TT,AT
        PPPQ=(l1[0]*l2[1])
        PQPP=(l1[1]*l2[0])
        PQQQ=(l1[1]*l2[2])
        QQPQ=(l1[2]*l2[1])

        #I2=(l1[0]*l2[0])+(l1[1]*l2[1])+(l1[2]*l2[2])    #AA,AA    AT,AT   TT,TT
        PPPP=(l1[0]*l2[0])
        PQPQ=(l1[1]*l2[1])
        QQQQ=(l1[2]*l2[2])

       ###Store these pairwise likelihoods in an array
#        PIBS[gg]=[I0,I1,I2]
        PIBS[gg]=[PPQQ,QQPP,PPPQ,PQPP,PQQQ,QQPQ,PPPP,PQPQ,QQQQ]  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT 


        ###The same as above, but this time we are also weighting the two genptype likelihoods by the probability of observing that genotype combination in the populatio
        ######Include probability of genotype
        #I0=(l1[0]*l2[2])+(l1[2]*l2[0])  #AA,TT    TT,AA
        PPQQg=l1[0]*l2[2]*(p**2)*(q**2)
        QQPPg=l1[2]*l2[0]*(q**2)*(p**2)

        #I1=(l1[0]*l2[1])+(l1[1]*l2[0])+(l1[1]*l2[2])+(l1[2]*l2[1])       #AA,AT    AT,AA   AT,TT   TT,AT
        PPPQg=l1[0]*l2[1]*(p**2)*(2*p*q)
        PQPPg=l1[1]*l2[0]*(2*p*q)*(p**2)
        PQQQg=l1[1]*l2[2]*(2*p*q)*(q**2)
        QQPQg=l1[2]*l2[1]*(q**2)*(2*p*q)

        #I2=(l1[0]*l2[0])+(l1[1]*l2[1])+(l1[2]*l2[2])    #AA,AA    AT,AT   TT,TT
        PPPPg=l1[0]*l2[0]*(p**4)
        PQPQg=l1[1]*l2[1]*(2*p*q)*(2*p*q)
        QQQQg=l1[2]*l2[2]*(q**2)*(q**2)

        
#        PIBS[gg]=[I0,I1,I2]
        PIBSg[gg]=[PPQQg,QQPPg,PPPQg,PQPPg,PQQQg,QQPQg,PPPPg,PQPQg,QQQQg]  #AATT,TTAA,AAAT,ATAA,ATTT,TTAT,AAAA,ATAT,TTTT 

    ###We transpose the matrices
    PIBSt=PIBS.transpose()
    PIBStg=PIBSg.transpose()

    ###identify the most likely genotype combination
    BestGT=PIBS.argmax(axis=1)

    BestIBS=np.zeros((nbSNPs,3),float)


    ###For each SNP, given the best genotype combination, pulls out the appropriate P(IBS|IBD) for all three IBS possibilities
    for gg in range(nbSNPs):
        BestIBS[gg]=IBS_all[BestGT[gg]][gg]

    ###ML optimization
    res=[]
    for gg in range(3):
        ok=0
        while ok==0:
            if gg==0:
                k1,k2=0.0,0.0
            else:
                k1,k2=random.random(),random.random()
            if k1+k2<=1.0:
                k0=1-(k1+k2)
                if 4*k2*k0<=k1**2:
                    k=np.array([k1,k2])
                    if kin(k,BestIBS,mask_mat)<>10E10:
                        ok=1
        temp=fmin(kin,k,args=(BestIBS,mask_mat),xtol=0.01,ftol=0.01,maxiter=None,maxfun=None,full_output=1, disp=0, retall=0, callback=None)                  
        res.append([temp[1],temp[0]])

    res2=[]
    for gg in range(3):
        ok=0
        while ok==0:
            if gg==0:
                k1,k2=0.0,0.0
            else:
                k1,k2=random.random(),random.random()
            if k1+k2<=1.0:
                k0=1-(k1+k2)
                if 4*k2*k0<=k1**2:
                    k=np.array([k1,k2])
                    if GLkin(k,PIBSt,IBS_all,mask_mat)<>10E10:
                        ok=1
        temp=fmin(GLkin,k,args=(PIBSt,IBS_all,mask_mat),xtol=0.01,ftol=0.01,maxiter=None,maxfun=None,full_output=1, disp=0, retall=0, callback=None)                  
        res2.append([temp[1],temp[0]])



    res.sort()
    res2.sort()
    
#    out= str(pw[g][0])+'\t'+str(pw[g][1])
    out=head[pw[g][0]]+'\t'+head[pw[g][1]]
    out=out+'\t'+str(round(1-(res[0][1][0]+res[0][1][1]),2))+'\t'+str(round(res[0][1][0],2))+'\t'+str(round(res[0][1][1],2))+'\t'+str(round(0.5*res[0][1][0]+res[0][1][1],3))
    out=out+'\t'+str(round(1-(res2[0][1][0]+res2[0][1][1]),2))+'\t'+str(round(res2[0][1][0],2))+'\t'+str(round(res2[0][1][1],2))+'\t'+str(round(0.5*res2[0][1][0]+res2[0][1][1],3))
    out=out+'\t'+str(len(mask_mat)-np.sum(mask_mat))+'\n'

    print out[:-1]

    file=open(filenameout,'a') 
    file.write(out)
    file.close()



