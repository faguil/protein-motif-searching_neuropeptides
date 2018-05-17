# Protein motif searching script for neuropeptide precursors

## Description

Perl script that scans one o more protein motif patterns in FASTA files (i.e., protein sequences). This scripts relies on PROSITE; therefore, it requires the external program pfscan. A binary executable file is given in this repository. This script has been tested on Linux Ubuntu 16.04, but it is probable that also works on others Linux operating systems.
Any further questions, please contact me at felipe.zoujiro@gmail.com

## Publication

This script is associated with the following publication:

- Thiel D, Franz-Wachtel M, Aguilera F & Hejnol A. (2018). Changes in the neuropeptide complement correlate with nervous system architecture in xenacoelomorphs. Molecular Biology and Evolution (Accepted)

## Before using protein_motif_searching.pl script

Install pfscan system-wide or in your $PATH

Copy the pfscan binary file to /usr/local/bin/

      sudo cp pfscan /usr/local/bin/

To run protein_motif_searching.pl script, a copy of the PROSITE database is needed. This database is also provided in this repository (i.e., prosite.dat). Keep the prosite.dat file in the same folder where the protein_motif_searching.pl script is located.

## Usage

      perl protein_motif_searching.pl -h

protein_motif_searching.pl [options] sequence-file(s)
protein_motif_searching version 1.1 options:
-h : this help screen

Input/Output:
  -p <string> : specify a protein motif pattern based on PROSITE patterns
  -o <string> : specify output format : gff

Other options (Use under your own risk, without warranty of any kind):
  -e <string> : specify the ID or AC (based on PROSITE database) of an entry in sequence-file
  --reverse   : randomize the sequence database by taking the reverse sequence of each individual entry

Note:
  * The sequence-file must be in FASTA format
  * There may be several -p arguments

## Examples

      perl protein_motif_searching.pl -p "G[KR](2)-x(2,35)-G[KR](2)-x(2,35)-G[KR](2)" o gff protein-file.fa > motif-pattern.txt

      perl protein_motif_searching.pl -p "G[KR](2)-x(2,35)-G[KR](2)-x(2,35)-G[KR](2)" -p "IIRIF" -o gff protein-file.fa > motif-pattern.txt
      

## Output

      perl protein_motif_searching.pl -p "G[KR](2)-x(2,35)-G[KR](2)-x(2,35)-G[KR](2)" -p "IIRIF" -o gff protein-file.fa > motif-pattern.txt

Gene.22058::Asco.Locus_10254.0_Transcript_1/0_Confidence_3_Length_333|FPKM|67.44::g.22058::m.22058	protein_motif_searching|v1.1	USER001	7	30	.	.	.	Sequence "GRKaikeqgs............................GKKtrkrnmkk...........................GKR" ; SequenceDescription "Gene.22058::Asco.Locus_10254.0_Transcript_1/0_Confidence_3_Length_333|FPKM|67.44::g.22058::m.22058"

Gene.16850::Asco.Locus_7367.0_Transcript_1/0_Confidence_2_Length_925|FPKM|0.51|BlastHit|swissprot|P20812|CP2A3_RAT::g.16850::m.16850	protein_motif_searching|v1.1	USER002	21	25	.	.	.	Sequence "IIRIF" ; SequenceDescription "Gene.16850::Asco.Locus_7367.0_Transcript_1/0_Confidence_2_Length_925|FPKM|0.51|BlastHit|swissprot|P20812|CP2A3_RAT::g.16850::m.16850"

Output explanation:

- The first column refers to sequence header from protein FASTA file containing the user-defined protein motif pattern (i.e., protein motif pattern used for searching).
- The second column refers to the name of the script and its version.
- The third column refers to arbitrary names for protein motif pattern used for searching (i.e., USER001 is the first protein motif pattern used for searching, while USER002 is the second motif pattern used for searching, etc.,).
- The fourth column refers to start position where user-defined protein motif pattern was found in the corresponding protein sequence.
- The fifth column refers to the end position where user-defined motif pattern was found in the corresponding protein sequence.
- The sixth, seventh and eighth columns are used for output format purposes. 
In fact, these columns are no informative.
- The ninth column refers to amino acid sequences (in CAPITAL letters) corresponding to the user-defined motif patterns. Depending on the user-defined motif pattern, you will obtain the amino acids (LOWER-CASE letter) surrounding the user-defined motif pattern used for searching.
- The tenth column refers to sequence description from the protein FASTA file containing the user-defined motif pattern (i.e., motif pattern used for searching).
