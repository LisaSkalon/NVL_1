# First homework
My 1st python homework. 2nd semester, Bioinformatics Institute.

## Needleman-Wunsh algorithm + argparse

***input:*** fasta file(string 1 and string 2), arguments  
***standart output:*** alignment + score matrix

### Usage and arguments
***usage***  
NVL_argparse.py [-h] -i Str [-s Str] [-g Int] [-p] [-n Str]  


***example***  
**input:** python NVL_argparse.py -i shortseqs.fasta -s 'My Matrix:' -g 5 -n 'output'  
**output:** file 'output.txt' with matrix:  
My Matrix:  
[[ 244.  234.  169.]  
 [ 234.  224.  159.]  
 [ 169.  159.   94.]]  


***optional arguments***  
  -h, --help    :          show this help message and exit  
  -i Str, --input Str   :  Input file (required) 
  -s Str, --string Str  :  String, printed before matrix  
  -g Int, --gap Int     :  Score for gap  
  -p, --print_flag     :   Print seq1 and seq 2 for me!  
  -n Str, --name Str   :   File name  
  
### Special features
Reading *fasta* file with *SeqIO*
```python
with open(in_file, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
```
### Authors
Lisa Skalon
