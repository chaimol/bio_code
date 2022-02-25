#!/usr/bin/bash
#作用是获取批量下载NCBI的数据脚本
case $1 in
	-h|--help)
		echo -e "
		Usage: \n
		bash getNCBI.bash -t SRR SRR12455818 SRR12455819 SRR12455820
		bash getNCBI.bash -t SRP SRP12345678
		bash getNCBI.bash -t ERP ERP12345678
		\n
		-t 指定输入的编号类型：支持SRR, ERR, SRP, ERP\n
		后面是对应的编号
		"
		exit 1
		;;
esac

#创建输出文件
random_str=`date +%s%N | md5sum|cut -c 1-5` #获取5位随机字符串
mkdir output_${random_str} #创建输出文件夹

NCBI_json="NCBI.${random_str}.json" #NCBI.XXXXX.json 文件的下载地址信息
info_list="info.${random_str}.list" #info.XXXXX.list json的解析结果
download_sh="download.${random_str}.sh" #download.XXXXX.sh 是ascp下载的命令
MD5="md5.${random_str}.status" #md5.XXXXX.status #MD5检测文件

#第1步，获取json文件
ffq -o $PWD/output_${random_str}/${NCBI_json} $@
if [ $? -eq 0 ];then
	echo "Output json info in $PWD/output_${random_str}/${NCBI_json}"
else
	echo "Error in step 1, get json file Error!"
	exit 1
fi
#第2步，解析json文件
json2tab.py $PWD/output_${random_str}/${NCBI_json} > $PWD/output_${random_str}/${info_list}
if [ $? -eq 0 ];then
	echo "Output tab info in $PWD/output_${random_str}/${info_list}"
else
	echo "Error in step 2, 解析 ${NCBI_json}.json 文件失败!请提交issue！"
	exit 1
fi
#第3步，获取下载地址
cat $PWD/output_${random_str}/${info_list}|awk '{print $3}'|sed 's/ftp:\/\/ftp.sra.ebi.ac.uk/ascp -QT -l 100m -P33001 -i \$HOME\/.aspera\/connect\/etc\/asperaweb_id_dsa.openssh era-fasp\@fasp.sra.ebi.ac.uk:/g;s/$/ ./g' >$PWD/output_${random_str}/${download_sh}
echo "Output download sh info in $PWD/output_${random_str}/${download_sh}"
#bash $PWD/output_${random_str}/${download_sh} & #直接开始下载

#第4步，检测下载文件的完整性
awk '{print $3,$5}' $PWD/output_${random_str}/${info_list}|rev|cut -d "/" -f1|rev|awk '{print $2,$1}' >$PWD/output_${random_str}/${MD5}
echo "Output md5 info in $PWD/output_${random_str}/${MD5}"
