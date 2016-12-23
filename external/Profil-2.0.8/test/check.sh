#! /bin/sh

echo "******************"
echo "*  Making check  *"
echo "******************"
p=0
w=0
f=0 
for t in $*; 
do 
	echo -n $t
	log=`echo $t | sed -e "s/\(.*\)\..*/\1.log/"`
	./$t >$log 2>&1; 
	r=$?
	if [ $r -eq 0 ]
	then 
		echo "     passed"
		p=`expr $p + 1`
	elif [ $r -eq 1 ]
	then 
		echo "     FAILED - please check $log"
		f=`expr $f + 1`
	else 
		echo "    passed with WARNING - please check $log"
		w=`expr $w + 1`
	fi
done; 

echo "************************************"
echo "  Passed $p of $#"
echo -n "  WARNING: $w "
if [ $w -gt 0 ]
then
	echo "- please check logs"
else
	echo
fi
echo -n "  FAILED : $f "
if [ $f -gt 0 ]
then
	echo "- please check logs"
else
	echo 
fi
echo "************************************"

