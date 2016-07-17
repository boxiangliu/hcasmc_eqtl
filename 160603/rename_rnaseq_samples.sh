#!/bin/bash
input=$1
dir=$2
while read line; do 
	line=($line)
	old=${line[0]}
	new=${line[1]}
	echo $old 
	echo $new 
	mv $dir/$old $dir/$new
done < $input