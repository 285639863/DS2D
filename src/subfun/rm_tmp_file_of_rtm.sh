#!/bin/sh
#Author: Duan Xinbiao,SGRI,2011.09.09
#Usage: ./rm_tmp_file_of_rtm.sh    Before using this script on one node, you must rsh to this node firstly.I wanted to execute the script to many nodes at same time, unfortunately failed.


    tmp='./dxbtmp'
  
    j=10
    ee=99
    while [ $j -le $ee ];do
        ls -lrt /scr02/dxb_rtm_ws*_isanp_$j* |awk '{ print $9} ' >>$tmp
        let j=j+1
        rm -f `cat $tmp`
        rm -f $tmp
    done 

    j=0
    ee=9
    while [ $j -le $ee ];do
        ls -lrt /scr02/dxb_rtm_ws*_isanp_$j* |awk '{ print $9} ' >>$tmp
        let j=j+1
        rm -f `cat $tmp`
        rm -f $tmp
    done 

