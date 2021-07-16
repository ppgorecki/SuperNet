
# Test on 5 species; 100 networks (50 pairs)
# embnet.py access needed

SPECIES=10
RET=4

SL=/tmp/s.log
EL=/tmp/e.log

for i in {1..100} 
do
	#nt1
	N=$( embnet.py -n "rand:$SPECIES:$RET" -pn )
	G=$( embnet.py -g "rand:$SPECIES:0" -pg )

	echo $G 
	echo $N

	#supnet -n "$N" -ed	
	echo "supnet -g \"$G\" -n \"$N\" -eo "
	S=$( supnet -g $G -n "$N" -eo -CDC | cut -f1 -d' '  )

	echo "==="
	E=$( embnet.py -n $N -g $G -pa -cc | grep -v "Net" | cut -f3 -d' ' | sort -n  | head -1 )

	echo $S $E

	if [[ $E != $S ]]
	then
		echo $N		
		echo $G
		echo "supnet -g \"$G\" -n \"$N\" -eo "
		break
	fi
done

rm -f $SL $EL