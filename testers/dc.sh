

# Test on 10 species; 100 trees (50 pairs)
# embnet.py access needed

supnet -r100 -A10 | xargs -n2 | while read ST GT
do
	echo "supnet -s \"$ST\" -g \"$GT\" -CDC -ec "
	A=$( supnet -s $ST -g $GT -CDC -ec )
	B=$( embnet.py  -g "$GT" -s"$ST" -pc -cc | cut -f3 -d' ' )
	
	if [[ $A != $B ]]
	then
		echo ERR in $ST $GT Supnet $A ref $B
		echo "supnet -s $ST -g $GT -CDC -ec"
		exit -1
		break
	else
		echo OK $ST $GT $A
	fi

done