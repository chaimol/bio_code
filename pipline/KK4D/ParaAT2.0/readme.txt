Readme of ParaAT (Oct. 4, 2014; version 2.0) 

1. Before running
Please make sure the following programs are put into the global directory and can be accessible from any working directory.
	(1) ParaAT.pl;
	(2) Epal2nal.pl;
	(3) Multiple sequence aligner (clustalw2 | t_coffee | mafft | muscle).

Sequence IDs in amino acid file, nucleotide file, and homologous group file should be the same. 

2. How to run
The test is on a small dataset containing 100 human-mouse orthologs.
	test.cds: coding sequences for human and mouse
	test.pep: protein sequences for human and mouse
	test.homologs: 100 pairs for human-mouse orthologs
	proc: number of processors, which is set at 6 and allowable for change during running

The command to run on this dataset is:
	ParaAT.pl -h test.homologs -n test.cds -a test.pep -p proc -o output

For help information and details parameter setting, please type:
	ParaAT.pl -help

3. Download & Update
Please visit <http://cbb.big.ac.cn> for downloads and updates. 

4. Contact information

Zhang Zhang, Ph.D.

CAS Key Laboratory of Genome Sciences and Information
Beijing Institute of Genomics, Chinese Academy of Sciences (CAS)
No.1 Beichen West Road
Chaoyang District, Beijing 100101
China
	
Email:	zhangzhang@big.ac.cn
Tel.:	+86(10)8409-7261
Fax:	+86(10)8409-7845
Web:	http://cbb.big.ac.cn
