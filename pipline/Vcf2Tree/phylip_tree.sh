#!/bin/bash
#目的：自动化生成phylip需要的par文件
#运行方法：bash phylip_tree.sh sample.phy sample_name

#sample_name是输出文件前缀
if [ $# -eq 0 ] || [ $# -eq 1 ];then
    echo "Usage:
        bash phylip_tree.sh sample.py sample_name"
        exit 1
fi

#定义输入文件
sample=$1  #phy文件
simple=$2  #输出结果文件前缀

#定义输出par函数
function make_par(){
#cat seqboot.par
echo "$sample
R
1000
Y
9" >$simple.seqboot.par
#cat dnadist.par
echo "$simple.seqboot.out
T
2.3628
M
D
1000
2
Y" >$simple.dnadist.par
#cat neighbor.par
echo "$simple.dnadist.out
M
1000
9
Y" >$simple.neighbor.par
# cat consense.par
echo "$simple.nei.tree
Y">$simple.consense.par
}


###par文件参数讲解
<<!
#cat seqboot.par
$sample #设定输入.phy文件的名称，否则输入默认的名为infile的文件
R #选择bootstrap
1000 #设置bootstrap的值，即重复的replicate的数目，通常使用1000或者100，注意此处设定好后，后续两步的M值也为1000或者100
Y #yes确认以上设定的参数
9 #设定随机参数，输入奇数值。

#cat dnadist.par
$simple.seqboot.out #本程序的输入文件
T #选择设定Transition/transversion的比值
2.3628 #比值大小
M #修改M值
D #修改M值
1000 #设定M值大小
2 #将软件运行情况显示出来
Y #确认以上设定的参数

#cat neighbor.par
$simple.dnadist.out #本程序的输入文件
M
1000  #设定M值大小
9 #设定随机数，输入奇数值
Y #确认以上设定的参数

# cat consense.par
$simple.nei.tree  #本程序的输入文件
Y #确认以上设定的参数
!

#定义生成tree文件的函数
function get_tree(){
seqboot<$simple.seqboot.par && mv outfile $simple.seqboot.out && \
dnadist<$simple.dnadist.par && mv outfile $simple.dnadist.out && \
neighbor<$simple.neighbor.par && mv outfile $simple.nei.out && mv outtree $simple.nei.tree  &&  \
consense<$simple.consense.par && mv outfile $simple.cons.out && mv outtree $simple.constree
}

#执行函数
make_par
get_tree
