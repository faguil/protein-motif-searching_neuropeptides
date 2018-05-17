# Protein motif searching script for neuropeptide precursors

## Description

Perl script that scans one o more protein motif patterns in FASTA files (i.e., protein sequences). This scripts relies on PROSITE; therefore, it requires the external program pfscan. A binary executable file is given in this repository. This script has been tested on Linux Ubuntu 16.04, but it is probable that also works on others Linux operating systems.

## Publication

This script is associated with the following publication:

- Thiel D, Franz-Wachtel M, Aguilera F & Hejnol A. (2018). Changes in the neuropeptide complement correlate with nervous system architecture in xenacoelomorphs. Molecular Biology and Evolution **(Accepted)**

## Before using protein_motif_searching.pl script

Install pfscan system-wide or in your $PATH

Copy the pfscan binary file to /usr/local/bin/

      sudo cp pfscan /usr/local/bin/

To run protein_motif_searching.pl script, a copy of the PROSITE database is needed. This database is also provided in this repository (i.e., prosite.dat.gz). First, you need to unzip the prosite.dat.gz 

      gunzip prosite.dat.gz

Keep the prosite.dat file in the same folder where the protein_motif_searching.pl script is located. 

## Usage/Examples

      perl protein_motif_searching.pl -p "G[KR](2)-x(2,35)-G[KR](2)-x(2,35)-G[KR](2)" o gff protein-file.fa > motif-pattern.txt

      perl protein_motif_searching.pl -p "G[KR](2)-x(2,35)-G[KR](2)-x(2,35)-G[KR](2)" -p "IIRIF" -o gff protein-file.fa > motif-pattern.txt
      

## Output

      perl protein_motif_searching.pl -p "G[KR](2)-x(2,35)-G[KR](2)-x(2,35)-G[KR](2)" -p "IIRIF" -o gff protein-file.fa > motif-pattern.txt

Gene.22058::Asco.Locus_10254.0_Transcript_1/0_Confidence_3_Length_333|FPKM|67.44::g.22058::m.22058	protein_motif_searching|v1.1	USER001	7	30	.	.	.	Sequence "GRKaikeqgs............................GKKtrkrnmkk...........................GKR" ; SequenceDescription "Gene.22058::Asco.Locus_10254.0_Transcript_1/0_Confidence_3_Length_333|FPKM|67.44::g.22058::m.22058"

Gene.16850::Asco.Locus_7367.0_Transcript_1/0_Confidence_2_Length_925|FPKM|0.51|BlastHit|swissprot|P20812|CP2A3_RAT::g.16850::m.16850	protein_motif_searching|v1.1	USER002	21	25	.	.	.	Sequence "IIRIF" ; SequenceDescription "Gene.16850::Asco.Locus_7367.0_Transcript_1/0_Confidence_2_Length_925|FPKM|0.51|BlastHit|swissprot|P20812|CP2A3_RAT::g.16850::m.16850"

Output explanation:

- The first column refers to sequence header from protein FASTA file containing the user-defined protein motif pattern (i.e., protein motif pattern used for searching).
- The second column refers to the name of the script and its version.
- The third column refers to arbitrary names for protein motif pattern used for searching (i.e., USER001 is the first protein motif pattern used for searching, while USER002 is the second protein motif pattern used for searching, and so on).
- The fourth column refers to start position where user-defined protein motif pattern was found in the corresponding protein sequence.
- The fifth column refers to the end position where user-defined protein motif pattern was found in the corresponding protein sequence.
- The sixth, seventh and eighth columns are used for output format purposes. **In fact, these columns are no informative.**
- The ninth column refers to amino acid sequences (in CAPITAL letters) corresponding to the user-defined protein motif patterns. Depending on the user-defined protein motif pattern, you will obtain the amino acids (LOWER-CASE letter) surrounding the user-defined protein motif pattern used for searching.
- The tenth column refers to sequence description from the protein FASTA file containing the user-defined protein motif pattern (i.e., motif pattern used for searching).

## Motif pattern matching nomenclature

Given that this script relies on PROSITE, its nomenclature is used for creating/generating the user-defined protein motif patterns used for searching.

Below critical information is given for creating/generating the user-defined protein motif patterns.

**Pattern syntax:**

The standard IUPAC one letter code for the amino acids is used in PROSITE.
The symbol 'x' is used for a position where any amino acid is accepted.
Ambiguities are indicated by listing the acceptable amino acids for a given position, between square brackets '[ ]'. For example: [ALT] stands for Ala or Leu or Thr.
Ambiguities are also indicated by listing between a pair of curly brackets '{ }' the amino acids that are not accepted at a given position. For example: {AM} stands for all any amino acid except Ala and Met.
Each element in a pattern is separated from its neighbor by a '-'.
Repetition of an element of the pattern can be indicated by following that element with a numerical value or, if it is a gap ('x'), by a numerical range between parentheses. 

**Examples:**

x(3) corresponds to x-x-x
x(2,4) corresponds to x-x or x-x-x or x-x-x-x
A(3) corresponds to A-A-A

When a pattern is restricted to either the N- or C-terminal of a sequence, that pattern respectively starts with a '<' symbol or ends with a '>' symbol. 

**Notes:**

Ranges can only be used with with 'x', for instance 'A(2,4)' is not a valid pattern element.

Ranges of 'x' are not accepted at the beginning or at the end of a pattern unless resticted/anchored to respectively the N- or C-terminal of a sequence, for instance 'P-x(2)-G-E-S-G(2)-[AS]-x(0,200)' is not accepted but 'P-x(2)-G-E-S-G(2)-[AS]-x(0,200)>' is.

If your pattern does not contain any ambiguous residues, you don't need to specify separation with '-'. 
Example: M-A-S-K-E can be written as MASKE. 

**Additional examples:**

[AC]-x-V-x(4)-{ED}	[Ala or Cys]-any-Val-any-any-any-any-{any but Glu or Asp}
<A-x-[ST](2)-x(0,1)-V	Ala-any-[Ser or Thr]-[Ser or Thr]-(any or none)-Val at the N-terminal of the sequence
<{C}*>	No Cys from the N-terminal to the C-terminal 
i.e. All sequences that do not contain any Cys.
IIRIFHLRNI	Ile-Ile-Arg-Ils-Phe-His-Leu-Arg-Asn-Ile

## License

Given this script relies on PROSITE, the Prosite.pm module is needed. This module is a Copyright (C) 2001-2006 to the Swiss Institute of Bioinformatics. It is released under the terms of the GNU General Public License, available at http://www.gnu.org/copyleft/gpl.html. This perl module is also provided in this repository and should be kept in the same folder where the protein_motif_searching.pl script is located.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Any further questions, please contact me at felipe.zoujiro@gmail.com
