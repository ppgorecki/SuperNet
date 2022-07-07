import pandas as pd
from subprocess import check_output, CalledProcessError, STDOUT
import os

Src="sim.csv"
cubedata = pd.read_csv(Src, sep=",")
env = dict(os.environ)

# make output dir: bbresults
check_output(f"mkdir -p bbresults", shell=True, universal_newlines=True, env=env)
makepdfs=""

# to make pdfs uncomment
# makepdfs="j"

f=open("sim_results.csv","w")
f.write("pair,n,r,network_nr,t,cost,time,minrt,naivecnt,naivetime,dpcnt,dptime\n")
for index,i in cubedata.iterrows():
	# take a gene_tree and a network

	gtr = i['gene_tree']
	net = i['network']
	network_nr = i['network_nr']
	n = i['n']
	r = int(i['r'])

	print(f"pair={index}\tn={n}\tr={r}")

	# tests for various -t thresholds
	for t in range(5,r+2):

		res = check_output(f"supnet -g '{gtr}' -n '{net}' -t{t} -ebk{makepdfs} -CDC", shell=True, universal_newlines=True, env=env)

		cost,time,minrt,naivecnt,naivetime,dpcnt,dptime=res.replace("\n","").split(" ")
		print(f"\tt={t}\tcost={cost}\ttime={time}\tminrt={minrt}\tnaivecnt={naivecnt}\tnaivetime={naivetime}\tdpcnt={dpcnt}\tdptime={dptime}")
		f.write(f"{index},{n},{r},{network_nr},{t},{cost},{time},{minrt},{naivecnt},{naivetime},{dpcnt},{dptime}\n")
						
		if makepdfs:
			check_output(f"dot -Tpdf bb.dot -o bbresults/{index}.t{t}.pdf ", shell=True, universal_newlines=True, env=env)			

f.close()


	




