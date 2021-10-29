
# https://doi.org/10.1016/j.jcss.2020.06.001
# Counts of one-component TCNs with k reticulations on [n] 1 ≤ k < n ≤ 8.
#kn	2	3	4   	 5     		 
#1	2	21	228 	 2805  		 
#2	-	42	1272	 30300 		 
#3	-	-	2544	 154500		 
#4	-	-	-		 309000		 

for n in 2 3 4 5; do
	for k in 1 2 3 4; do
		if [[ $k -lt $n ]]
		then
			echo $n $k
			supnet -r1000000 -R$k -A$n -eu >/dev/null 2>$n.$k.cnt &
		fi
	done
done
wait 

for i in *cnt; do echo $i $( cat $i) ; done



