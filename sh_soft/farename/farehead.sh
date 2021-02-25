#!/bin/bash
#farehead 主要是对fa文件头部提取字符。
#Usage: farehead.sh demo2.fa
  #用于修改fa文件的头部。注意修改里面的参数

case $1 in
	-h|--help)
echo -e "Usage \n
purpose:for exchange fasta file head. \n
   -h help;\n
   -eh input.fa \n
	
Example
  
  # use a new head file
  ./farehead.sh -eh demo2.fa >demo2.new.fa
 
 
Note:
可能根据你的fa的头部的实际情况，调整第37行的处理方式
"
		;;
		-eh | --exchange_head)
			exchange_head $2
		  ;;
		*)
	  echo "Input info is wrong,please test input -h or --help for the help info!"
	  ;;
esac	  

exchange_head(){
  cat $1|while read line;
  do
   echo $line|grep "=" >cache.txt
  if [ $? -eq 0 ];                  #判断函数运行返回值，等于0，则成功，不等于0，则检查用户的输入，告知错误原因！
  then
          echo $line|cut -d "=" -f4|cut -d " " -f1|sed 's/^/\>/'  #需要调整的行
 else
          echo $line
 fi
 done
 rm -rf cache.txt
}