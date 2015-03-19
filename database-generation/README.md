# database-generation

## Howto: use this script

1) Download the .vcf.gz file from the [Exome Aggregation Consortium](http://exac.broadinstitute.org/) and decompress it. Only 0.2 and 0.3 releases have been tested. The file will take around 20Go afterwards. 

2) Import the database-generation project into an IDE. Here we will assume you will use Eclipse.

3) In the Main.java file, set the number of threads the script can use. Set also the number of samples you want to generate, and if you want one big file with all the samples or 1 file by sample.

4) In the same file, also indicate the path to the .vcf previously downloaded, and the directory where you want to generate the new .vcf files.

5) In the panel "Run Configurations" for the project, go to the tab "Arguments" and set the additional VM arguments: "-d64 -Xms9g -Xmx9g". Set the memory according to your hardware obviously, it will also work with 3Go of Ram or any value (but you should allocate at least 2Go).

6) That's all, just run the project now.