
C=DC
if [[ $1 ]]; then C=$1; fi
LASTSP=j

echo Testing $C

while true
do

	ST=$( gsevol.py -g"randbin(a-$LASTSP)" -egP )
	GT=$( gsevol.py -g"randbin(a-$LASTSP,a-$LASTSP,a-$LASTSP)" -egP )
	#GT=$( gsevol.py -g"randbin(a-$LASTSP)" -egP )
	GT=$( gsevol.py -g"randbin(a-$LASTSP,a-$LASTSP,a-j)" -egP )
	# echo ST=\"$ST\"
	# echo GT=\"$GT\"
	# echo "supnet -s \$ST -g \$GT -C$C -ec"  
	S=$( supnet -s $ST -g $GT -C$C -ec )

	if [[ $C = D ]]
	then	
		G=$( gsevol.py -g "$GT" -s "$ST" -eGmih0 -w'(1,0,0)' -C'dltcostformat="!c"' )
	elif [[ $C = DC ]]
		then		
		G=$( embnet.py  -g "$GT" -s"$ST" -pc -cc | cut -f3 -d' ' )
	elif [[ $C = L ]]
		then		
		G=$( gsevol.py -g "$GT" -s "$ST" -eGmih0 -w'(0,1,0)' -C'dltcostformat="!c"' )		
	elif [[ $C = DL ]]
		then		
		G=$( gsevol.py -g "$GT" -s "$ST" -eGmih0 -w'(1,1,0)' -C'dltcostformat="!c"' )		
	else
		echo "??" 
		break
	fi


	if [[ $S != $G ]]
	then
		echo ERR in $ST $GT Supnet $S ref $G
		break
	else
		echo OK $ST $GT $S
	fi

	#break

done