#!/bin/bash -x
#
i=$1


#replace this date with start date
year=$((`date --date="9 Apr 2011 +$i days" +%-Y`))
mon=$((`date --date="9 Apr 2011 +$i days" +%-m`))
day=$((`date --date="9 Apr 2011 +$i days" +%-d`))
today=$((`date --date="9 Apr 2011 +$i days" +%Y%m%d`))
var=$(( $i + 1 ))
tomorrow=$((`date --date="9 Apr 2011 +$var days" +%Y%m%d`))
nyear=$((`date --date="9 Apr 2011 +$var days" +%-Y`))
var=$(( $i - 1 ))
yesterday=$((`date --date="9 Apr 2011 +$var days" +%Y%m%d`))
pyear=$((`date --date="9 Apr 2011 +$var days" +%-Y`))


echo $yesterday
echo $today

#radars that you want to search through
for rad in "fhw" "gbr"
do

	############################################################################
	#this block of code makes sure that 3 consecutive days of filtered
	#files exist centered on the target date (the one you set above)
	#you will have to change this code, but that is what you need to make it do
	############################################################################
	
  if test $i -eq 0
  then
		if [ ! -e "/data/fit/$yesterday.${rad}.fitexbfc" ] || [ ! -s "/data/fit/$yesterday.${rad}.fitexbfc" ]
		then
			cp /sd-data/$pyear/rawacf/$rad/$yesterday* .
			bunzip2 *.bz2
			gunzip -f *.gz
			cat $yesterday*$rad*rawacf > $yesterday.${rad}.rawacf
			make_fitex2 -new $yesterday.${rad}.rawacf > /data/fit/$yesterday.${rad}.fitex
			fitexfilter /data/fit/$yesterday.${rad}.fitex > /data/fit/$yesterday.${rad}.fitexbfc
			rm $yesterday*
		fi
  fi

  if test $i -lt 1
  then
  	if [ ! -e "/data/fit/$today.${rad}.fitexbfc" ] || [ ! -s "/data/fit/$today.${rad}.fitexbfc" ]
		then
			cp /sd-data/$year/rawacf/$rad/$today* .
			bunzip2 *.bz2
			gunzip -f *.gz
			cat $today*$rad*rawacf > $today.${rad}.rawacf
			make_fitex2 -new $today.${rad}.rawacf > /data/fit/$today.${rad}.fitex
	#     mv $today.${rad}.rawacf /data/raw
			fitexfilter /data/fit/$today.${rad}.fitex > /data/fit/$today.${rad}.fitexbfc
			rm $today*
		fi
  fi

	if [ ! -e "/data/fit/$tomorrow.${rad}.fitexbfc" ] || [ ! -s "/data/fit/$tomorrow.${rad}.fitexbfc" ]
	then
		cp /sd-data/$nyear/rawacf/$rad/$tomorrow* .
		bunzip2 *.bz2
		gunzip -f *.gz
		cat $tomorrow*$rad*rawacf > $tomorrow.${rad}.rawacf
		make_fitex2 -new $tomorrow.${rad}.rawacf > /data/fit/$tomorrow.${rad}.fitex
	#   mv $tomorrow.${rad}.rawacf /data/raw
		fitexfilter /data/fit/$tomorrow.${rad}.fitex > /data/fit/$tomorrow.${rad}.fitexbfc
		rm $tomorrow*
	fi
	#########################################################################################
	#end the section of code that makes sure the 3 days of files exist
	#########################################################################################

	
	#concatenate the 3 files into a single file
  cat /data/fit/$yesterday.${rad}.fitexbfc /data/fit/$today.${rad}.fitexbfc /data/fit/$tomorrow.${rad}.fitexbfc > 2_day.${rad}.fitexbfc

  #execute the search, give info on the target date (year, mon, day)
  ./dopsearch -new -year $year -day $day -mon $mon /data/fit/$today.${rad}.fitexbfc > $today.$rad.gscat


	#check that the file is not nonzero size
  FILESIZE=$(stat -c%s "$today.$rad.gscat")
  echo $FILESIZE
  if test $FILESIZE -gt 0
  then
		#move the file to a new directory
    mv $today.$rad.pscat /data/raw/pscat/$today.$rad.pscat
  fi

	#clean up files
  rm $today*

done
