#!/bin/bash

## usage: 
## script1.sh [file_list]

i=0

echo "Files Processed:"

for file in $(ls -1 *transcript_idPA.txt)
##for file in $(cat $1)
do
	
	if [ $i == 0 ] 
	then
		##cat $file > raw_counts.counts
		awk '{print $1"\t"$2}' $file > raw_counts.counts
	else 
		awk '{print $2}' $file > temp.xls
		paste raw_counts.counts temp.xls > comb_copy.xls
		cp comb_copy.xls raw_counts.counts
	fi

	echo "$file"
	((i++))

done

rm temp.xls comb_copy.xls
