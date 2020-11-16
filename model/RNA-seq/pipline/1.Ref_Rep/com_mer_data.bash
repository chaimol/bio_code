#!/bin/bash
##运行方式
#bash com_mer_data.bash wangll.trans.csv 
##程序目的：提取转录组的转录本的产出的文件的基因名

##直接运行程序的时候指定文件名
filename1=$1

###判断循环，如果运行时候直接有参数，则直接后续运行，若无参数 ，则要求用户输入文件名。
if [ -n "$1" ]; then
    echo "包含第1个参数，继续运行。"
else
    ##读取用户输入的文件名，并储存在filename变量中
	echo "请输入文件名"
	read -p "input the transcript  count filename:" filename1   
fi

cat $filename|awk -F "," '{print $1}' | awk -F "_" '{print $1}' >gene_id
cat $filename|awk -F "," '{print $1}' >trans_id
cat $filename|tr "," " " >new_trans
##获取第二列到最后一列的值
#cat $filename|awk -F "," '{for (i=2;i<NF;i++)printf("%s ",$i);print ""}'|tr "," " " >count_info
##paste是连接多个文件 -d 指定连接符   tr是转换字符，第一个值是文中的字符，第二个是需要更改为的字符。
paste -d " " gene_id new_trans|tr " " "," >new_trans_count.csv
rm -rf gene_id trans_id new_trans
echo "complete exchange id and add gene_id!"


<<!
#说明

#####程序运行方法#######
###bash com_mer_data.bash file1.csv file2 1 1 csv
##程序参数说明
##参数1和参数2是两个文件名，是必须的。
##参数3和4是选填，指定需要比对的列，默认是第一列。
##参数5是指定文件的分割符号的类型，默认是csv格式。

###程序用于合并两个具有相同列的文件
###
###

filename2=$2
ref1=$3
ref2=$4

if [ -n "$2" ]; then
    echo "包含第2个参数，继续运行。"
else
    ##读取用户输入的文件名，并储存在filename变量中
	echo "请输入文件名"
	read -p "input the transcript  count filename:" filename2   
fi




show_usage="args: [-f1 , -f2 , -r1 ,-r2,-t, -h]\
                                  [--file1=, --file2=, --row1=, --row2=,--type,--help]"
#参数
# 本地仓库目录
opt_localrepo=""

# git仓库url
opt_url=""

# 备份目录
opt_backupdir=""

# web目录
opt_webdir=""

GETOPT_ARGS=`getopt -o l:r:b:w: -al local-repository:,repository-url:,backup-dir:,webdir: -- "$@"`
eval set -- "$GETOPT_ARGS"
#获取参数

while [ -n "$1" ]
do
        case "$1" in
                -l|--local-repository) opt_localrepo=$2; shift 2;;
                -r|--repository-url) opt_url=$2; shift 2;;
                -b|--backup-dir) opt_backupdir=$2; shift 2;;
                -w|--webdir) opt_webdir=$2; shift 2;;
                --) break ;;
                *) echo $1,$2,$show_usage; break ;;
        esac
done

if [[ -z $opt_localrepo || -z $opt_url || -z $opt_backupdir || -z $opt_webdir ]]; then
        echo $show_usage
        echo "opt_localrepo: $opt_localrepo , opt_url: $opt_url , opt_backupdir: $opt_backupdir , opt_webdir: $opt_webdir"
        exit 0
fi




<<
#!/bin/sh
#说明
show_usage="args: [-l , -r , -b , -w]\
                                  [--local-repository=, --repository-url=, --backup-dir=, --webdir=]"
#参数
# 本地仓库目录
opt_localrepo=""

# git仓库url
opt_url=""

# 备份目录
opt_backupdir=""

# web目录
opt_webdir=""

GETOPT_ARGS=`getopt -o l:r:b:w: -al local-repository:,repository-url:,backup-dir:,webdir: -- "$@"`
eval set -- "$GETOPT_ARGS"
#获取参数
while [ -n "$1" ]
do
        case "$1" in
                -l|--local-repository) opt_localrepo=$2; shift 2;;
                -r|--repository-url) opt_url=$2; shift 2;;
                -b|--backup-dir) opt_backupdir=$2; shift 2;;
                -w|--webdir) opt_webdir=$2; shift 2;;
                --) break ;;
                *) echo $1,$2,$show_usage; break ;;
        esac
done

if [[ -z $opt_localrepo || -z $opt_url || -z $opt_backupdir || -z $opt_webdir ]]; then
        echo $show_usage
        echo "opt_localrepo: $opt_localrepo , opt_url: $opt_url , opt_backupdir: $opt_backupdir , opt_webdir: $opt_webdir"
        exit 0
fi

!
