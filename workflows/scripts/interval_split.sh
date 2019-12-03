#!/bin/bash -x

IN_FILE=$1 ;
O_BASE=$2 ;

#make sure tmpdir exists for sort
mkdir -pv $TMPDIR ;

EXTENSION=`echo $IN_FILE|grep -Po '\.[^\.]*$'|tr -d "."` ;
echo "EXTENSION FOUND TO BE $EXTENSION" ;

#loop over chromosomes
for CHROM in `cut -f1 $IN_FILE |grep -Pv '^@'| sort|uniq`; do

	#define output
	SPLIT_FILE=$O_BASE.$CHROM.$EXTENSION ;
	echo "To write to $SPLIT_FILE"  ;

	#write headers if they exist
	cat $IN_FILE|grep -P '^@' > $SPLIT_FILE ;

	#write the intervals, subset with awk
	cat $IN_FILE | awk -F "\t" -v CHROM_COMP="$CHROM" '{if($1==CHROM_COMP) print $0}' >> $SPLIT_FILE ;

	#see if output has any intervals!
	NUM_INTS=`cat $SPLIT_FILE|awk -F "\t" -v CHROM_COMP="$CHROM" '{if($1==CHROM_COMP) print $0}' |cut -f2,3|grep -Pc '\d+\s+\d+'` ;
	if [ "$NUM_INTS" -gt "0" ] ; then
		# a valid interval file
		echo "At least one mutect interval detected on sequence $CHROM " ;
		#save integers to file for subsequent iteration for indexing
	else
		# no intervals detected!
		echo "No mutect intervals detected on sequence $CHROM .  So deleting data split on the sequence . " ;
		rm -v $SPLIT_FILE ;
	fi ;

done ;