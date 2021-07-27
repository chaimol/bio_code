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
- Step3: `KK4D.sh -h` for help

# Use
#### for coline analysis
`KK4D.sh coline`

#### for kaks calculator
`KK4D.sh kaks`

#### for 4DTv calculator
`KK4D.sh 4DTV`

#### from gff3 get bed 
`KK4D.sh bed`

#### from gff3 and protein.fa ,get the protein sequence
`KK4D.sh pep`

#### from gff3 and cds.fa, get the longest transcript sequence
`KK4D.sh cds`

#### from gff3 cds.fa  protein.fa ,get 1 or 2 species all the above information.
`KK4D.sh all `

# Warning
Please make sure that the input sequence IDs of cds.fa and protein.fa are the same.
Such as:
`head Ath.protein.fa`
```
>AT3G05780.1 pep chromosome:TAIR10:3:1714941:1719608:-1 gene:AT3G05780 transcript:AT3G05780.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:LON3 description:Lon protease homolog 3, mitochondrial [Source:UniProtKB/Swiss-Prot;Acc:Q9M9L8]
MMPKRFNTSGFDTTLRLPSYYGFLHLTQSLTLNSRVFYGARHVTPPAIRIGSNPVQSLLL
FRAPTQLTGWNRSSRDLLGRRVSFSDRSDGVDLLSSSPILSTNPNLDDSLTVIALPLPHK
PLIPGFYMPIHVKDPKVLAALQESTRQQSPYVGAFLLKDCASTDSSSRSETEDNVVEKFK
VKGKPKKKRRKELLNRIHQVGTLAQISSIQGEQVILVGRRRLIIEEMVSEDPLTVRVDHL
```
`head Ath.cds.fa`
```
>AT3G05780.1 cds chromosome:TAIR10:3:1714941:1719608:-1 gene:AT3G05780 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:LON3 description:Lon protease homolog 3, mitochondrial [Source:UniProtKB/Swiss-Prot;Acc:Q9M9L8]
ATGATGCCTAAACGGTTTAACACGAGTGGCTTTGACACGACTCTTCGTCTACCTTCGTAC
TACGGTTTCTTGCATCTTACACAAAGCTTAACCCTAAATTCCCGCGTTTTCTACGGTGCC
CGCCATGTGACTCCTCCGGCTATTCGGATCGGATCCAATCCGGTTCAGAGTCTACTACTC
```
`head Ath.gff3`
```
1       TAIR10  chromosome      1       30427671        .       .       .       ID=chromosome:1;Alias=CP002684.1,Chr1,NC_003070.9
1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;Name=NAC001;biotype=protein_coding;description=NAC domain-containing protein 1 [Source:UniProtKB/Swiss-Prot%3BAcc:Q0WV96];gene_id=AT1G01010;logic_name=araport11
1       araport11       mRNA    3631    5899    .       +       .       ID=transcript:AT1G01010.1;Parent=gene:AT1G01010;biotype=protein_coding;transcript_id=AT1G01010.1
1       araport11       five_prime_UTR  3631    3759    .       +       .       Parent=transcript:AT1G01010.1
1       araport11       exon    3631    3913    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon1;constitutive=1;ensembl_end_phase=1;ensembl_phase=-1;exon_id=AT1G01010.1.exon1;rank=1
1       araport11       CDS     3760    3913    .       +       0       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
1       araport11       exon    3996    4276    .       +       .       Parent=transcript:AT1G01010.1;Name=AT1G01010.1.exon2;constitutive=1;ensembl_end_phase=0;ensembl_phase=1;exon_id=AT1G01010.1.exon2;rank=2
1       araport11       CDS     3996    4276    .       +       2       ID=CDS:AT1G01010.1;Parent=transcript:AT1G01010.1;protein_id=AT1G01010.1
```
In this example,the cds and protein IDs is "AT3G05780.1",so we must be set the "key" and "type" in config.ini. "key" should be set is "transcript_id" and "type" should be set is "mRNA".
If this example,if the "key" is set to "ID", then the output IDs of Ath.bed will be "transcript:AT1G01010.1" is different with the protein.fa and cds.fa IDs "AT1G01010.1".So this will course Error.

Before each script is run, it will check whether the output file of the previous step exists. So if the output bed format is incorrect, you can manually adjust it and still rename it to the name of the program output. Then run subsequent commands. If the latter step fails, be sure to delete the failed file.
# Output info

 
# Update information
#### 2021.3.18 release the Version 0.01
#### 2021.3.19 update the Version to 0.02
Update info:
1. The V0.01 has to much bug.
2. All commands in this version have been tested and run normally.
3. This version does not include a visualization module, and the next version may add the visualization module.




