
# Test on 5 species; 100 networks (50 pairs)
# embnet.py access needed

SPECIES=10
RET=9

SL=/tmp/s.log
EL=/tmp/e.log

for i in {1..100} 
do
	#tree child test		
	#N=$( embnet.py -n "rand:$SPECIES:$RET" -pn )

	#nt1
	N=$( embnet.py -n "nt1:$SPECIES:$RET" -pn )

	#generic network - now with errors, also in embnet.py
	#N="(a,(c,((b)#B,((#B)#A,#A))))"
	#N=$( embnet.py -n "gen:$SPECIES:$RET" -pn )

	#echo $N	

	#supnet -n "$N" -ed	
	supnet -n "$N" -ed | supnet -A $SPECIES -S- -eS  | sort | uniq > $SL
	
	#embnet.py -n $N -pd
	embnet.py -n $N -pd | grep -v "Network" | supnet -A $SPECIES -S- -eS  | sort | uniq > $EL
	
	echo $N  $( wc -l $SL ) $( wc -l $EL )

	if ! diff $SL $EL
	then
		echo DIFF $N
		cat $SL
		echo "---"
		cat $EL
		break
	fi
done

rm -f $SL $EL