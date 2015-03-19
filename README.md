# cgs-benchmarks
Description of the benchmarking studies peformed in this CGS system. This module is part of the [**CGS**](https://github.com/jpoullet2000/cgs) project. The goal of this module is to describe different benchmarks that have been carried out within the CGS project to assess the efficiency of the system and its scalability. 

## Tools

### database-generation
Super-optimized script using the files provided by the [Exome Aggregation Consortium](http://exac.broadinstitute.org/) to generate vcf files from the main vcf file available [here](http://exac.broadinstitute.org/downloads). The script allows to create one vcf by sample or one vcf for all samples. The user can decide how many samples he wants to be generated. Code will obviously be faster if you decide to put all samples in one vcf file, as the main bottleneck is the speed of writing data to your disk.