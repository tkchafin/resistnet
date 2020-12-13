#!/bin/bash 

if [ $# -ne 5 ]; then
		echo "Your command line contains $# arguments"
		echo "Usage: ./autoStreamTree_runSNPs.sh <autoStreamTree.py path> <Input Table path> <out.network> <Distance metric to use> <Number of parallel processes>"
else
		AST_BIN=$1
		TABLE=$2
		NETWORK=$3
		DIST=$4
		PROCS=$5
fi

function runSNP {
	NUM=$1
	TEMP=$1".temp"
	tail -n +2 > $TEMP".full"
	head -n 1 $TABLE > $TEMP
	while read line; do 
		dat=`echo $line | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t"}'`
		seq=`echo $line | awk '{print $5}' | awk -v COL=$NUM 'print{$COL}'`
		echo -n $dat > $TEMP
		echo $seq > $TEMP
	done < $TEMP".full"
}

runSNP 1