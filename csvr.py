#!/usr/bin/python3

import pandas as pd
from subprocess import check_output, CalledProcessError, STDOUT
import os
import sys
from pathlib import Path

force=False
dotopt=""

def bashu(comm):
	return check_output(comm,shell=True, universal_newlines=True, env=env)	

def bashuclean(comm):
	res=bashu(comm)
	# ignore debug's
	for l in res.split("\n"):
		if l and l[0]=='#':
			print(l)

	res = [ l for l in res.split("\n") if l and l[0]!='#' ]
	return res[0]

# Usage:
# input_file.csv minrt|minlb|bfs [limitexpr]

printcommand = False

if len(sys.argv)<3:
	print(f"""
Usage:
{sys.argv[0]} [pdf|png] [printcommand] [force] [initscore PRC% ] input_file.csv [minrt|minlb|bfs]+ [limitexpr]

pdf|png - make pdf or png files with BB tree
limitexpr - any python expr with n, r, t and index variables
printcommand - print supnet command
force - force computation
initscore PRC% - test BB with init best score: best_score*(PRC/100)+1; PRC=100 - optimal_scoretest (optcut)

Examples:
{sys.argv[0]} sim.csv bfs 
{sys.argv[0]} sim.csv bfs "n==40"
{sys.argv[0]} rnd.csv minrt "n==40 and r==20 and t==0"
{sys.argv[0]} rnd.csv png minrt "t>5 and index==20"
{sys.argv[0]} initscore 100 sim.csv minlb 

{sys.argv[0]} sim.csv minlb minrt bfs
for i in 100 120 140 160 180 200; {sys.argv[0]} initscore $i rnd.csv minlb minrt bfs "(t>3) and (t%2==0)"; end
{sys.argv[0]} rnd.csv minlb minrt bfs
{sys.argv[0]} initscore 100 rnd.csv minlb minrt bfs
		""")
	exit(0)

tp=""
nlim = None
args = sys.argv[1:]

env = dict(os.environ)
initscoretest = False		

while args and args[0] in ('pdf','png','printcommand','force','initscore'):

	com = args.pop(0)

	if com in ('pdf','png'):
		dotopt="j"	
		picouputdir="bbresults"
		bashu(f"mkdir -p {picouputdir}")
		pic=com
		continue

	if com == "printcommand":
		printcommand = True		
		continue
		scoreteststep=int(args.pop(0))
		scoretestcnt=int(args.pop(0))
		continue

	if com == "force":
		force = True		
		continue

	if com == "initscore":
		initscoretest = True		
		initscoretestprc=int(args.pop(0))
		continue

	break

import socket
hostname=socket.gethostname()

Dest=f"multitests_{hostname}"
Src=args.pop(0).rsplit(".")[0]
pdatadir=f"{Dest}/{Src}/pdata"

tps=[]
while args and args[0] in ('minrt','minlb','bfs'):
	tps.append(args.pop(0))

if args: nlim=args.pop(0)

cubedata = pd.read_csv(Src+".csv", sep=",")
	
def runtp(tp):

	if tp=="minrt":
		outputfile=f"minrt"		
	elif tp=="minlb":
		outputfile=f"minlb"
	elif tp=="bfs":	
		outputfile=f"bfs"
	else:
		print("Expected arg: bfs, minrt or minlb")
		exit(-1)

	if initscoretest: 
		outputfile=f"{outputfile}_initscore_{initscoretestprc}"	
		
	supnet=f"supnet_{tp}"
	print (outputfile,supnet)

	if not Path(supnet).is_file():
			print(bashu(f"make {supnet}"))

	print(bashu(f"mkdir -p {pdatadir}"))
		
	fmain=open(f"{Dest}/{Src}/{outputfile}.csv","w")
	fmain.write("pair,n,r,network_nr,t,cost,time,minrt,naivecnt,naivetime,dpcnt,dptime\n")

	errfiles = []

	cnt=0
	tmpfiles = []
	for index,i in cubedata.iterrows():
		# take a gene_tree and a network
		gtr = i['gene_tree']
		net = i['network']
		network_nr = i['network_nr']
		n = int(i['n'])
		r = int(i['r'])

		tvals = range(0,r+2)
		if nlim is not None:		
			try:   
				lst=[]
				for t in tvals:

					if eval(nlim, globals(), locals()): lst.append(t)                         
				tvals=lst			

			except NameError as e:                    
				print(e)
				exit(-1)

		if not tvals: continue

		cnt+=1

		csvf = f"{pdatadir}/{index}.{outputfile}.csv"
		
		tmpfl = f"/tmp/{index}.{Src}.{outputfile}.tmp"
	
		if not force and Path(csvf).is_file():
			print(f"#{cnt} FILE {csvf} exists")

		elif not force and Path(tmpfl).is_file():		

			print(f"#{cnt} FILE {tmpfl} exists")

			tmpfiles.append(csvf)

			continue

		else:

			f=open(tmpfl,"w")
			
			# print(f"#{cnt} pair={index}\tn={n}\tr={r}")
			
			# compute init score 			
			if initscoretest:
				# here
				comm=f"supnet -g '{gtr}' -n '{net}' -ebk -CDC"			
				res = bashuclean(comm)
				bestscore = int(res.split(" ")[1])
				#print(res)
				initscoredce = int(bestscore*initscoretestprc/100)
				print(f"initscoretest {initscoretestprc}% with cost={bestscore} -> initcost={initscoredce}")
				

			# tests for various -t thresholds
			for t in tvals:

				comm=f"{supnet} -g '{gtr}' -n '{net}' -t{t} -ebk{dotopt} -CDC"
				if initscoretest: comm+=f" -b{initscoredce}"
				res = bashuclean(comm).split(" ")

				if printcommand: print(comm)

				cost,costdce,time,minrt,naivecnt,naivetime,dpcnt,dptime=res
				# print(f"#{index}.{tp}\tt={t}\tcost={cost}\tcostdce={costdce}\ttime={time}\tminrt={minrt}\tnaivecnt={naivecnt}\tnaivetime={naivetime}\tdpcnt={dpcnt}\tdptime={dptime}")
				f.write(f"{index},{n},{r},{network_nr},{t},{cost},{time},{minrt},{naivecnt},{naivetime},{dpcnt},{dptime}\n")
								
				if dotopt:
					bashu(f"dot -T{pic} bb.dot -o {picouputdir}/{index}.t{t}.{pic}")			
					print(f"{picouputdir}/{index}.t{t}.{pic} created")

			f.close()

			print(bashu(f"mv {tmpfl} {csvf}"))

			print(f"{csvf} created")
		
		# check correctness
		f = open(csvf,"r")
		l = f.read()
		err = 0
		for i in l.split():
			if len(i.split(","))!=12: 
				errfiles.append(csvf)
				err = 1
		if not err: 
			fmain.write(l)
		f.close()

	if tmpfiles:
		print("Tmp - files:")
		for i in tmpfiles:
			print(i)

	if errfiles:	
		print("Err - files:")
		for i in errfiles:
			print(i)	

	fmain.close()
	print(outputfile+".csv created")


for tp in tps:
	runtp(tp)


