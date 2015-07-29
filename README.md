lcMLkin (Maximum Likelihood Estimation of Relatedness)
--------------------------------------------

How to Compile:

*Note:  This code requires a C++11 compliant compiler.  g++ >= 4.7.3 and clang >= 3.4 should be OK*

1.) `> git clone https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation.git`

2.) `> cd maximum-likelihood-relatedness-estimation`

3.) `> make`

Now, you should have an executable in this directory called `lcmlkin`.  You can run it with the `-h` option
to get help.  A typical run would look something like:

4.) `> lcmlkin -i input.vcf -o output.relate -g all -t 8`

This will run `lcmlkin` on the file `input.vcf`, estimate background allele frequencies assuming all individuals are unrelated, and using the variant of the algorithm that sums over all genotypes given their likelihoods (`-g all`).  It will write the output to the file `output.relate`, and will make use of up to 8 worker threads (`-t 8`).
