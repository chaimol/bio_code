# getNCBI.bash 用于获取下载NCBI的数据的地址和MD5值的文件

Author Email: `chaimol@163.com`

# 依赖软件：
- python3 模块 json,pandas,sys
```
pip install json
pip install pandas
```
- ffq [ffq github](https://github.com/pachterlab/ffq)

安装ffq 
  ```
  pip install ffq
  ```
- aspera的安装，请一定要安装这种方法安装，以确保最终的密钥地址是`$HOME\/.aspera\/connect\/etc\/asperaweb_id_dsa.openssh`
```
wget https://d3gcli72yxqn2z.cloudfront.net/connect_latest/v4/bin/ibm-aspera-connect_4.1.1.73_linux.tar.gz
tar -zxvf ibm-aspera-connect_4.1.1.73_linux.tar.gz
./ibm-aspera-connect_4.1.1.73_linux.sh
echo "export PATH=$HOME/.aspera/connect/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

# 安装本软件的方法
```
chmod 757 getNCBI.bash
chmod 757 json2tab.py
echo "export PATH=$PWD:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

# 使用方法
```
bash getNCBI -t SRR SRR12455818 SRR12455819 SRR12455820
bash getNCBI -t SRP SRP12345678
bash getNCBI -t ERP ERP12345678
```
目前测试过的包括SRR,ERR,ERP,SRP,这4种输入类型，其他的可能会有报错！

输入的参数包括：
- -t指定输入的编号类型：支持SRR, ERR, DRR, SRP, ERP, DRP, GSE, DOI
- 后续是编号信息

上面两个输入的参数和ffq的一致


# 输出文件（XXXXX是5位随机字符）
所有文件输出在output_XXXXX文件夹内
- NCBI.XXXXX.json #文件的下载地址信息
- info.XXXXX.list #输出的文件
- download.XXXXX.sh #是ascp下载的命令
- md5.XXXXX.status #MD5检测文件

# 注意事项
默认使用的是100m带宽的下载速度，如果你的带宽更大或更小，请修改
例如：
- 修改为30m `sed -i 's/-l 100m/-l 30m/g' download.*.sh` 
- 修改为1000m `sed -i 's/-l 100m/-l 1000m/g' download.*.sh`

# 后续分析
- 直接使用`bash download.*.sh &`即可开始下载数据
- 数据下载完成后，`md5sum -c md5.*.status`即可检测下载文件的完整性








