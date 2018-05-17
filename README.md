# Protein motif searching script for neuropeptide precursors

## Description

Perl script that scans one protein motif patterns in protein FASTA files. This scripts relies on PROSITE; therefore, it requires the external program pfscan. A binary executable file is given in this repository. This script has been tested on Linux Ubuntu 16.04, but it is probable that also works on others Linux operating systems.

## Before using protein_motif_searching.pl script

Install pfscan system-wide or in your $PATH

Copy the pfscan binary file to /usr/local/bin/

      sudo cp pfscan /usr/local/bin/

