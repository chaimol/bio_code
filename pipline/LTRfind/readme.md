# LTRfind: a pipeline for get LTR and LAI ( base on [LTR_retriever](https://github.com/oushujun/LTR_retriever))
# 1. Install require software by conda
```
conda create -n LTR_retriever
conda activate LTR_retriever
conda install -y -c conda-forge perl perl-text-soundex
conda install -y -c bioconda cd-hit repeatmasker
git clone https://github.com/oushujun/LTR_retriever.git
./LTR_retriever/LTR_retriever -h
conda install -c bioconda seqkit
chmod 757 LTRfind
./LTRfind -h
```
you also can add  the path of `LTRfind` to the `~/.bashrc`.

# 2.how to use it
It is mainly used for the identification of diploid and polyploid LTR and the calculation of LAI value.
### Purpose:
Get the LTR_RTs in the genome and calculate the LAI value of the genome at the same time. 
for Diploid use (`-D|Diploid`)
for Polyploid (up to octoploid) use (`-P|Polyploid`)

### Usage: 
```
	-h|--help 
	-D|Diploid Species genome.fa 
	-D|Diploid Species genome.fa ChrString 
	-P|Polyploid -abbr Species -g genome.fa -chr ChrString -p1 Ghir_A -p2 Ghir_D 
	-P|Polyploid -abbr Species -g genome.fa -p1 Ghir_A -p2 Ghir_D 
	-V | --version
```
### Example:
```
	LTRfind -D A.thaliana /path/A.thaliana.fa 
	LTRfind Diploid A.thaliana /path/A.thaliana.fa chr
	LTRfind -P -abbr Species -g genome.fa -chr ChrString -p1 Ghir_A -p2 Ghir_D 
	LTRfind -P -abbr Species -g genome.fa -chr ChrString -p1 Ghir_A -p2 Ghir_D -p3 Ghir_C 
	LTRfind -P -abbr Species -g genome.fa -chr ChrString -p1 Ghir_A -p2 Ghir_D -p3 Ghir_C -p4 Ghir_B
	LTRfind Polyploid -abbr Species -g genome.fa -p1 Ghir_A -p2 Ghir_D
```
### Note:
In run `-D|Diploid` ，if you input the ChrString , will extract the Chromosome sequence for LTR and LAI (not analysis the scaffold sequence).
Input genome file must be `.fa`,not support `.fa.gz`.
**the genome.fa sequence ID should be less then 15 character. if your file is report over 15 character, can use `seqkit seq -i genome.fa >new.genome.fa` to extract the ID and delete the info after the ID**
if your input `genome.fa` `ChrID` may be like `>1`,in this case, you cannot use the `ChrStr` field to extract chromosome data and filter scaffold. 
	It is recommended to use `LTRfind` after filtering by yourself. The ID length of genome.fa cannot exceed 15 characters. 
# 3. Input
### Diploid: 
- speciesname 
- genome.fa 
- ChrStr (can be omitted)
### Polyploid:
- speciesname 
- genome.fa 
- ChrStr (If you do not provide this parameter, it is assumed that the genome file you provide only contains the sequence of the chromosome.)
- p1
- p2
- p3
- p4
## Example:
#### Polyploid:
for `Ghirsutum.fasta`
`head Ghirsutum.fasta`
```
>Ghir_A01
TAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAA
ACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACC
CTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTA
AACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAC
CCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAACCCTAAACCCTAAACCC
TAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAA
ACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACC
CTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTA
AACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAC
...
>Ghir_D01
AAACCCTAAACCCTAAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAACCCTAAACCCTAAACCCTAAACCCTA
AACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAAC
CCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCT
AAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAA
CCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCCTAAACCC
```
should be run like this:
`LTRfind -P -abbr G.hirsutum -g Ghirsutum.fasta -chr Ghir_ -p1 Ghir_A -p2 Ghir_D`
#### Diploid
`head Arabidopsis_thaliana.fa`
```
>chr01 dna:chromosome chromosome:TAIR10:1:1:30427671:1 REF
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAAT
CTTTAAATCCTACATCCATGAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTT
CTCTGGTTGAAAATCATTGTGTATATAATGATAATTTTATCGTTTTTATGTAATTGCTTA
TTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCTTGT
GGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAA
```
`bash LTRfind -D A.thaliana Arabidopsis_thaliana.fa`or `bash LTRfind -D A.thaliana Arabidopsis_thaliana.fa chr`
# 4. Output
All the output is the same with [LTR_retriever](https://github.com/oushujun/LTR_retriever).
For Chinese, there are detail info [links](https://www.jianshu.com/p/ed289822c825) in Chinese.
# 5. Author & Version
Build date:2021.06.17
Last  update: 2021.06.18
Version: 0.02
Author: Mol Chai
Email: chaimol@163.com
