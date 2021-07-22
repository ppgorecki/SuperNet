
# Test on 5 species; 100 networks (50 pairs)
# embnet.py access needed

SPECIES=60
RET=8

for i in {1..1} 
do
	#nt1
	N=$( embnet.py -n "rand:$SPECIES:$RET" -pn )
	
	echo $N
	supnet -n "$N" -o -eM

	# echo "==="
	# E=$( embnet.py -n $N -g $G -pa -cc | grep -v "Net" | cut -f3 -d' ' | sort -n  | head -1 )

	# echo $S $E

	# if [[ $E != $S ]]
	# then
	# 	echo $N		
	# 	echo $G
	# 	echo "supnet -g \"$G\" -n \"$N\" -eo "
	# 	break
	# fi
done
