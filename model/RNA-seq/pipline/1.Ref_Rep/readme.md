this file is for analysis RNA-seq,wangll have 3 repeat .The species studied is maize (**Zeamays**)

有生物学重复，研究物种是玉米。

- wangll.s  是运行主程序的引导程序，主要是运行于PBS/qsub 任务管理服务器集群上的任务提交程序。运行方法是`qsub wangll.s`,然后会在后台自动运行主程序wangll.bash.
- wangll.bash 是实际运行的主程序，是从下机数据质控到最终产出表达矩阵的脚本。
- rename.bash 是在下机数据的存放目录运行，主要是备份数据和重命名下机数据的文件名。
- com_mer_data.bash 是用于提取转录组产出的文件的矩阵的基因列的基因名称