#请按照下面顺序安装对应的软件
#1.add path
chmod 757 calculate_4DTV_correction.pl
chmod 757 axt2one-line.py
chmod 757 run.sh
chmod 757 getKaKs.sh

echo "export PATH=${PWD}:\$PATH" >>~/.bashrc


#2.use conda install requirement software
conda create -n mmdetection python=3.7
conda activate mmdetection
conda install -y jcvi seqkit mafft

echo "export PATH=${PWD}/KaKs_Calculator2.0/bin/Linux:\$PATH" >>~/.bashrc
echo "export PATH=${PWD}/ParaAT2.0:\$PATH" >>~/.bashrc
source ~/.bashrc
