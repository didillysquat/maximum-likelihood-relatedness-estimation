Readme for SNPbam2vcf.py
========================

*  This python script that takes a list of bamfiles and performs variant calling such that the output is in the format required for lcMLkin to calculate pairwise relatedness.

*  Please note it is very slow, and is only designed for calling genotypes for specific SNPs. I strongly do not recommend using this as a general variant caller.

*  We are working on a better parser (this was an experimental one developed to test lcMLkin), but for now, if using 100 samples with 100K loci, it may need to run overnight.

*  [Pysam](http://pysam.readthedocs.org/en/latest/) and numpy must be installed in the version of python used

*  Bam files must be indexed. The SNP file must have the tab separated fields in the following order: chromosome, position (one-based co-ordinates), reference allele, alternate allele

*  Written (poorly) by Krishna Veeramah

General usage is: 
```
./SNPbam2vcf.py <bamlist> <fileoutname> <target_SNP_file>
```


Using the example files provided (data is in `data/SNPbam2vcf_data/`), it would be run using the following command line:
```
./SNPbam2vcf.py bam.list out_play.vcf SNP.list 
```
