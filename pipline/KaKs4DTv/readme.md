# Readme.md
### this is a pipline for analysis of coline genes,KaKs and 4DTv .
#### It can analysis 2 species or 1 species.

# Input 
###  Prepare Input file(can be normal or *.gz)
genome.gff3 
genome.pep.fa 
genome.cds.fa

# Require
conda,linux

# Install
- Step1:Install software
```
bash Install.sh
source ~/.bashrc
```
- Step2: Modify the configuration file `config.ini`
- Step3: `run.sh` for help

# Use
#### for coline analysis
`run.sh coline`

#### for kaks calculator
`run.sh kaks`

#### for 4DTv calculator
`run.sh 4DTV`

#### from gff3 get bed 
`run.sh bed`

#### from gff3 and protein.fa ,get the protein sequence
`run.sh pep`

#### from gff3 and cds.fa, get the longest transcript sequence
`run.sh cds`

#### from gff3 cds.fa  protein.fa ,get 1 or 2 species all the above information.
`run.sh all `



# Update information
2021.3.18 release the Version 0.01

