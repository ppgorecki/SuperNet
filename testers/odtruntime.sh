

GENETREES=30000
NETSIZE=100
GENESIZE=26 #max 26
RET=8

N=$( embnet.py -n "rand:$NETSIZE:$RET" -pn )

time supnet -r$GENETREES -A$GENESIZE | supnet -G- -n $N -eo > /dev/null