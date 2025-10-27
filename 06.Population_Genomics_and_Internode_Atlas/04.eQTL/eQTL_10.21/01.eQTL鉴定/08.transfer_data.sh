#!/bin/bash

# Created on 2025-03-27
# @author: xiazhongqiang92@163.com

# 从曙光传到小服务器
# 进入小服务器
/data0/agis_xiazhongqiang/Biosoft/fasttrans/rayfile-c -a harbin02.hpccube.com -P 65012 -u agis_xiazq -w 9380ec0f4bd549fa41-9100-4c5c-acb0-528e7b33b8d2 -no-meta -symbolic-links follow -retry 10 -retrytimeout 30 -o download -s '/project/02.YZ/R4.PopVar/pop_expression/04.eQTL/07.pre.all.data' -d /data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/01.data


