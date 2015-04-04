# Highlander-benchmarks

## Introduction
This script has the goal to establish a connection with the highlander database to be able to insert vcf files and execute sql queries.

## Howto: First connection with OpenVPN

1. Download openvpn-install-2.3.6-I001-x86_64.exe [here](https://openvpn.net/index.php/open-source/downloads.html)

2. Decompress the keys and certificates previously received.

3. Copy the files to "C:\Program Files\OpenVPN\config\"

4. Execute the OpenVPN as administrator. Then, right-lick on the OpenVPN GUI icon and select "connect".

## Howto: Insert vcf into Highlander

1. Download the vcf_import.py script and move it next to dbBuilder.jar (this jar comes with Highlander). Download also db-generation.jar.

2. Establish a connection with OpenVPN

3. Give some information about the ram, disk, etc. the vcf_import.py can use (edit the first lines in the python script)

4. If you have dbBuilder\_fast2.jar (around 100x or even 1000x times faster than dbBuilder.jar), install the database you received with it. If you have a SSD disk it will be around 1000x times faster, a hard-drive will only give you a speed up of 100x. If you use dbBuilder\_fast2.jar and a hard-drive, a file of 150Mo and 55k variants will be analyzed in ~3-5 minutes instead of 8h with the current dbBuilder.jar. 

5. Launch vcf_import.py, it will automatically download the public database, create random vcf from it and upload them to highlander