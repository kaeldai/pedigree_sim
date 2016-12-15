# Pedigree simulator

Simulates sequencing data for a given pedigree - including germline mutations, somatic mutations, and library read errors. Takes in a pedigree (ped) file and a reference (fasta), and produces data by the following procedure:  

1. Generates pedigree founder gametic DNA based on the reference.  

2. Creates gametic DNA for the second generation, using meiosis transition of founders DNA plus the given mutation rate. Repeats until all members have been processed.
3. Creates somatic DNA by simulating mitotic transition from germline, with error rates.
4. Simulates a library read call with specified depth from somatic DNA, using a mixed Multinomial Dirichlet distribution.

The output produces the following files:
* A tad file containing all library reads at every site.
* A vcf file containing the true gametic and somatic genotype for every site.
* A csv file with various statistics.

## Installation
Currently on runs on posix, built in C++ and requires the following libraries:
* [htslib](https://github.com/samtools/htslib)
* [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [boost](http://www.boost.org/)
* [GNU Scientific Library](https://www.gnu.org/software/gsl/)

## Usage

```
./simulator [options] --ped <family.ped> --region <chr:begin-end> <reference.fasta>
```
* Options  
 * --read-depths [unsigned int] : number of library calls to make at each site (default 30)
 * --prefix [string] : prefix file name for output .tad, .vcf and .csv (uses family id in pedigree by default)
 * --theta [float] : population diversity for founders (default 0.001)
 * --ref-weight [float] : bias towards the reference at each site (default 1.0, no bias)
 * --nuc-freqs ["float,float,float,float"] : population frequency of A,C,G,T respectivily (default "0.3,0.2,0.2,0.2")
 * --mu [float] : germline mutation rate (default 1e-8)
 * --mu-somatic [float] : somatic mutation rate (default 1e-8)
 * --mu-library [float] : library preparation error rate (default 1e-8)
 * --gamma [...] : two sets of parameters used in the dirichlet multinomial models.
 * --germline-mutation [X:i,...] : Forces a gametic mutation for indivual X at location i in contig (starting at 0). For multiple gametic mutations use a comma separated string, ie "NA1:0,NA1:1,NA2:1000". Currently does not work with pedigree founders.
 * --somatic-mutation [X:i,...] : Forces a somatic mutation for indivual X at location i in contig (starting at 0). For multiple somatic mutations use a comma separated string, ie "NA1:0,NA1:1,NA2:1000".
 * --seed [unsigned int] : Specify a seed for the random generation of germline and somatic mutations. If 0 (default) then seed will be based on the time function."
