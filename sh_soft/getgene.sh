
#主要是对vcf和gtf进行操作，提取基因。
#提取所有基因

function get_all_gene(gtf,output){
	gtf=$1
	output=$2
	cat $gtf|awk '$3=="gene" {print $0}' >$output
}

function get_vcf_region(input_vcf,region,output_vcf){
	input_vcf=$1
	region=$2
	output_vcf=$3
	
	if expr "-3" : '-\?[0-9]\+$' >/dev/null; then
		echo "染色体"; 
	elif expr $region ;then
		echo "正常区域"
	fi
	
	
	test=10:51549498-51590734
	chr=cat $region|cut -d ':' -f1
	left_region=cat $region|awk -F ":" '{print $2}'|cut -d '-' -f1
	right_region=cat $region|cut -d '-' -f2
	cat $input_vcf| awk '$1 == "$chr" && $2 >= $left_region && $2 <= $right_region {print $0} ' >$output_vcf
} 

function get_gtf_gene(gtf,region,output){
	gtf=$1
	region=$2
	output=$3
	
	test=10:51549498-51590734
	chr=cat $region|cut -d ':' -f1
	left_region=cat $region|awk -F ":" '{print $2}'|cut -d '-' -f1
	right_region=cat $region|cut -d '-' -f2
	cat $gtf| awk '$3=="gene" && $1 == "$chr" && $4 >= $left_region && $5 <= $right_region {print $0} ' >$output
}

case $1 in
	-h|--help)
echo -e "Usage \n
  -h help;\n
  get_all_gene inputgtf outputgtf \n
  get_vcf_region inputvcf region outputvcf (region example: chr10:1549498-1590734 or 10:1549498-1590734 or chr10 or 10) \n
  get_gtf_gene inputgtf region output (region example: chr10:1549498-1590734 or 10:1549498-1590734 or chr10 or 10)
	
Example
  
  #  get all gene from gtf.  
  getgene get_all_gene hg19.gtf hg19.all_gene.gtf
  
  #  get some region gene from vcf
  getgene get_vcf_region RIL334.pass.vcf 10:51549498-51590734
  
  # get chr10 vcf from genome vcf 
  getgene get_vcf_region RIL334.pass.vcf chr10
  
  #get some region gene from gtf
  getgene get_gtf_gene hg19.gtf 10:51549498-51590734 chr10.candidate.gtf
  
  
  
Note:
Region can be chr10 or 10 or chr10:1549498-1590734 or 10:1549498-1590734 ,but the  chromosome must be same with the input gtf/vcf file.   
"
		;;
	get_all_gene)
		get_all_gene $2 $3
		  if [ $? -eq 0 ];                  #判断函数运行返回值，等于0，则成功，不等于0，则检查用户的输入，告知错误原因！
		  then
			echo "success!"
		  elif ! [ -n "$3" ];
		  then
			echo "please add the output filename!"
		  else
			echo "Please check your input info ,exchange file not success"
		  fi
		;;
				*)
	  echo "input info is wrong,please test input -h find the help info!"
	  ;;
	get_gtf_gene)
		get_gtf_gene $2 $3 $4
		  if [ $? -eq 0 ];                  #判断函数运行返回值，等于0，则成功，不等于0，则检查用户的输入，告知错误原因！
		  then
			echo "success!"
		  elif ! [ -n "$4" ];
		  then
			echo "please add the output filename!"
		  else
			echo "Please check your input info ,exchange file not success"
		  fi
		;;
	get_vcf_region)
		get_vcf_region $2 $3 $4
		  if [ $? -eq 0 ];                  #判断函数运行返回值，等于0，则成功，不等于0，则检查用户的输入，告知错误原因！
		  then
			echo "success!"
		  elif ! [ -n "$4" ];
		  then
			echo "please add the output filename!"
		  else
			echo "Please check your input info ,exchange file not success"
		  fi
		;;
		*)
	  echo "Input info is wrong,please test input -h or --help for the help info!"
	  ;;
esac	  
