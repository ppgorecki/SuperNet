
# embnet.py access needed

SPECIES=10  # species 
RET=8      # reticulations 
GENETREES=4 # number of gene trees 

# OPTIONS
#
# run for NNI:
# testers/hillclimb.sh N 
#
# 


#Generate gene trees:
G=$( embnet.py -g "rand:$SPECIES:0" -pg )
for i in $( seq $GENETREES ); do
	G="$G;"$( embnet.py -g "rand:$SPECIES:0" -pg )
done

echo $G 

for i in {1..10} 
do
	#nt1
	N=$( embnet.py -n "rand:$SPECIES:$RET" -pn )  # start from random
	
	echo $N

	#supnet -n "$N" -ed	
	echo "supnet -g \"$G\" -n \"$N\" -on$1 "
	time supnet -g $G -n "$N" -on$1

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
