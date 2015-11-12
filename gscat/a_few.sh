#!/bin/bash -x
#number of days you want to run the search on
for (( i = 0 ; i < 1; i++ ))
do
	./compare_search.sh $i
done
