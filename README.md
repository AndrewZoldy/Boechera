# Boechera
Boechera project


### Length parser

Presented script allows you to obtain lengths of contigs ( total and by each species ) aligned on anythig instead plants.

#### Launch parameters

Current version has 3 launch parameter:

  * **-i (--input_file)** - path to input .fasta file with contigs
  * **-a (--alignment)** - path to input file with alignments (cutted)
  * **-a (--output_file)** - path to output file (.txt formatted)
  
#### Results

As result you will obtain an .txt file with the following structure: 
>*Genus Species* : *summary length of contigs aligned on current Sp.*
 
 *Genus2 Species2* : *summary length of contigs aligned on current Sp.*
 
 ...
 
 unknown :  *summary length of contigs which did not align on anything*
 
 Total : *summary length of each contigs which passed through plants-filter*
 
 Total without unknowns : *summary length of each contigs which passed through plants-filter without unknown contigs*
