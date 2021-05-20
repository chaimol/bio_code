#!/bin/bash
chmod 757 Vcf2Tree
chmod 757 calc_4dTv_in_eff_vcf.py
chmod 757 phylip_tree.sh
chmod 757 vcf2phylip.py
chmod -R 757 ${PWD}/fastme-2.1.5/binaries
chmod -R 757 ${PWD}/iqtree-1.6.12/bin
chmod -R 757 ${PWD}/phylip-3.697/exe
chmod -R 757 ${PWD}/VCF2Dis-1.30/bin
#add software path to ~/.bashrc
echo "PATH=${PWD}/fastme-2.1.5/binaries:\$PATH" >>~/.bashrc
echo "PATH=${PWD}/iqtree-1.6.12/bin:\$PATH" >>~/.bashrc
echo "PATH=${PWD}/phylip-3.697/exe:\$PATH" >>~/.bashrc
echo "PATH=${PWD}/VCF2Dis-1.30/bin:\$PATH" >>~/.bashrc

echo "PATH=${PWD}:\$PATH" >>~/.bashrc
source ~/.bashrc
