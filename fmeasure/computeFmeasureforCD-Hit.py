#!/usr/bin/env python
import sys
import os

fastafilename = sys.argv[1] # the fasta file
namepos = int(sys.argv[2]) # the position of the feature in the header of the sequences used for comparisons
threshold = float(sys.argv[3])
endthreshold = 0
if len(sys.argv) > 4:
	endthreshold = float(sys.argv[4])
step = 0.01
if len(sys.argv) > 5:
	step = float(sys.argv[5])

def GetClasses(fastafilename):
	classes = []
	speciesnames = []
	#read the fasta file
	fastafile= open(fastafilename)
	for line in fastafile:
		if line.startswith('>'):
			words = line.split('|')
			seqid = int(words[0].strip('>'))
			speciesname = words[namepos]
			if speciesname in speciesnames:
				i=speciesnames.index(speciesname)
				refclass = classes[i] 
				refclass.append(seqid)
			else:
				speciesnames.append(speciesname)
				refclass=[]
				refclass.append(seqid)
				classes.append(refclass)
	return classes

def GetClusters(resultfilename):
	clusters = []
	#read the clustering result file
	cluster = []
	i=0
	resultfile = open(resultfilename)
	for line in resultfile:
		if line.startswith('>'):
			cluster = []
			clusters.append(cluster)
		else:
			words = line.split('>')
			seqid = int(words[1].split('|')[0])
			cluster.append(seqid)
	return clusters

def ComputeFmeasure(classes,clusters):
	#compute F-measure
	f=0
	n=0
	for refclass in classes:
		m = 0
		for cluster in clusters:
			i = len(set(refclass) & set(cluster))
			v = 2*i/(len(refclass) + len(cluster))
			if m < v:
				m=v
		n = n + len(refclass)
		f = f +	(len(refclass)*m)
	f = f /n 
	return f

if endthreshold < threshold:
	#cluster and compute F-measure
	command = "cd-hit-est -i " + fastafilename + " -n 10 -o result -c " + str(threshold)
	#run command
	os.system(command)
	classes = GetClasses(fastafilename)
	clusters = GetClusters("result.clstr")
	f = ComputeFmeasure(classes,clusters)
	print("F-measure: " + str(f))
else:
	resultfile = open("opt.txt","a")
	classes = GetClasses(fastafilename)
	#predict optimal threshold
	t = threshold
	while t <= endthreshold:
		#cluster and compute F-measure
		command = "cd-hit-est -i " + fastafilename + " -n 10 -o result -c " + str(t)
		#run command
		os.system(command)
		clusters = GetClusters("result.clstr")
		f = ComputeFmeasure(classes,clusters)
		resultfile.write(str(t) + '\t' + str(f) + '\n')
		t = t + step
	resultfile.close()
	os.system("rm result*")
		
