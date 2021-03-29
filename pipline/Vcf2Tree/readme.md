# Vcf2Tree 
Vcf2Tree is a pipline for make a evolutionary tree. It is the downstream analysis of GTAK4. 
Vcf2Tree is a tool to build a phylogenetic tree using SNP.vcf（or 4DTV.vcf）!

# Input file
- **sample.vcf** (The length of the material name grouped in vcf must be <= 10 characters. Characters > 10 will be automatically truncated to 10 characters.)
(if run 4DTV,the input.vcf must from snpEff!)
input.vcf must vcf,vcf.gz is not support!
- **output_prefix**  
- genome.fa (only required when run `4DTv`!)
# output file
-viq ${output}.min4.phy.contree and ${ouput}.min4.phy.treefile
-vdp ${output}.cons.vdp.tree
-vp ${output}.nei.tree and ${output}.contree
-vdf ${output}.vdf.tree
4DTv ${output}.4dtv.vcf

# Usage:
```
		#get tree by iqtree (run fast，can auto choose the best suitable Model) Highly recommended!!!
		Vcf2Tree -viq input.vcf output_prefix
		Vcf2Tree --vcf2iqtree input.vcf output_prefix
		
		#get tree by VCF2Dis and phylip (third fast)
		Vcf2Tree -vdp input.vcf output_prefix
		Vcf2Tree --vcf2Dis_phylip input.vcf output_prefix
		
		#get tree by VCF2Dis and fastme (second fast,Inconsistent with the results constructed by other methods )
		Vcf2Tree -vdf input.vcf output_prefix
		Vcf2Tree --vcf2Dis_fastme input.vcf output_prefix
		
		#get tree by phylip (slow！but accurate)
		Vcf2Tree -vp input.vcf output_prefix
		Vcf2Tree --vcf2phyliptree input.vcf output_prefix
		
		# get 4DTV.vcf
		Vcf2Tree 4DTv input.vcf output_prefix genome.fa
		
		# get 4DTV.vcf and tree 
		Vcf2Tree 4DTv input.vcf output_prefix genome.fa -viq
		Vcf2Tree 4DTv input.vcf output_prefix genome.fa -vdp
		Vcf2Tree 4DTv input.vcf output_prefix genome.fa -vp
		Vcf2Tree 4DTv input.vcf output_prefix genome.fa -vdf
		
		demo:
		Vcf2Tree -viq ${PWD}/demo/50k.vcf 50k
		Vcf2Tree -vdp ${PWD}/demo/50k.vcf 50k
		Vcf2Tree -vp ${PWD}/demo/50k.vcf 50k
		
		Notes:
		The output tree file can be visual by :
		Figtree https://github.com/rambaut/figtree/releases (base Java)
		ITOL https://itol.embl.de/upload.cgi(online)
		ggtree https://github.com/YuLab-SMU/ggtree (R packages)
```
# Require software
The require software have been put in this packages file.Just run `bash Install.sh`
# Install
```
bash Install.sh
source ~/.bashrc
```


