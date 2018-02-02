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
	clusterindexes = []
	i=0
	resultfile = open(resultfilename)
	for line in resultfile:
		words = line.split('\t')
		if len(words) > 1:
			clusterindex = int(words[1])
			seqid = int(words[8].split('|')[0])
			if clusterindex in clusterindexes:
				i = clusterindexes.index(clusterindex)
				clusters[i].append(seqid)
			else:
				cluster = []
				cluster.append(seqid)
				clusters.append(cluster)
				clusterindexes.append(clusterindex)
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

#sort the sequences
path = "/home/dvu/tools/uclust/"
command = path + "./uclust --sort " + fastafilename  + " --output seqs_sorted.fasta"
os.system(command)

if endthreshold < threshold:
	#cluster
	command = path + "./uclust --input seqs_sorted.fasta --uc result.uc --id " + str(threshold)
	#run command
	os.system(command)
	classes = GetClasses(fastafilename)
	clusters = GetClusters("result.uc")
	print(len(clusters))
	f = ComputeFmeasure(classes,clusters)
	print("F-measure: " + str(f))
else:
	resultfile = open("uclustopt.txt","a")
	classes = GetClasses(fastafilename)
	#predict optimal threshold
	t = threshold
	while t <= endthreshold:
		#cluster and compute F-measure
		#run command
		command = path + "./uclust --input seqs_sorted.fasta --uc result.uc --id " + str(t)
		os.system(command)
		clusters = GetClusters("result.uc")
		f = ComputeFmeasure(classes,clusters)
		resultfile.write(str(t) + '\t' + str(f) + '\n')
		t = t + step
	resultfile.close()
	os.system("rm result*")
		
