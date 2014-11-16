#!/bin/sh
# @Author: Marc
# @Date:   2014-11-15 23:48:11
# @Last Modified by:   Marc
# @Last Modified time: 2014-11-16 00:59:48

if [ $# -eq 0 ] || [ -z "$1" ] || [ "$#" -ne 3 ]
	then
	echo "Enter the path directory where are the songs to cut"
	echo "$0 <path> <start_time> <end_time>"
	exit 1
fi

if [ ! -d "$1" ]
	then
	echo "The directory does not exist"
	exit 2
fi

path=$1
files="$1/*.wav"
start_time=$2
end_time=$3

for file in $files
do
	if [[ -f $file ]]
		then
		filename=$(basename "$file")
		filename="${filename%.*}"
		echo "${filename}.wav is being trimmed..."
		./sox "$file" "${path}/${filename}_trimmed.wav" "trim" "$2" "$3"
		echo "done"
	fi
done