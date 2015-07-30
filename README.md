lcMLkin (Maximum Likelihood Estimation of Relatedness)
======================================================

[![Join the chat at https://gitter.im/COMBINE-lab/maximum-likelihood-relatedness-estimation](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/COMBINE-lab/maximum-likelihood-relatedness-estimation?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

How to Compile:

*Note:  This code requires a C++11 compliant compiler.  g++ >= 4.7.3 and clang >= 3.4 should be OK*

1.) `> git clone https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation.git`

2.) `> cd maximum-likelihood-relatedness-estimation`

3.) `> make`

Now, you should have an executable in this directory called `lcmlkin`.  You can run it with the `-h` option
to get help.  A typical run would look something like:

4.) `> lcmlkin -i input.vcf -o output.relate -g all -t 8`

This will run `lcmlkin` on the file `input.vcf`, estimate background allele frequencies assuming all individuals are unrelated, and using the variant of the algorithm that sums over all genotypes given their likelihoods (`-g all`).  It will write the output to the file `output.relate`, and will make use of up to 8 worker threads (`-t 8`).

Producing VCF input from BAM files
-----------------------------------

A preliminary implementation of a script to generate lcMLkin-ready VCF files from a set of indexed BAM files is available under the `src_python/SNPbam2vcf` directory of the repository ([here](https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation/tree/master/src_python/SNPbam2vcf)).  There, you will find more information on the script including usage examples.

