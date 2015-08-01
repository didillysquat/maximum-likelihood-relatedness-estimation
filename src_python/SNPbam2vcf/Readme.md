Readme for SNPbam2vcf.py
========================

*  This python script takes a list of bamfiles and performs variant calling such that the output is in the format required for lcMLkin to calculate pairwise relatedness.

*  Please note it is slow, and is only designed for calling genotypes for specific SNPs in smallish numbers (tens to hundreds) of individuals. I do not recommend using this as a general variant caller.

*  We are working on a better version (this was an experimental one developed to test lcMLkin), but for now, if using 100 samples with 100K loci with 20x mean coverage, this takes about 4 hours on one processor. Lower coverage datasets, which is what lcMLkin is designed for, will be quicker though, maybe just over an hour.

*  [Pysam](http://pysam.readthedocs.org/en/latest/) and [Numpy](http://www.numpy.org/) must be installed in the version of python used.

*  Bam files must be indexed. 

*  The SNP file must have the tab separated fields in the following order: chromosome, position (one-based co-ordinates), reference allele, alternate allele.

*  If you want to change things like mapping and base quality threshold, edit the python code under the section "Input arguments".
 
*  Written (poorly) by Krishna Veeramah

General usage is: 
```
./SNPbam2vcf.py <bamlist> <fileoutname> <target_SNP_file>
```


Using the example files provided (data is in `data/SNPbam2vcf_data/`), it would be run using the following command line:
```
./SNPbam2vcf.py bam.list out_play.vcf SNP.list 
```
