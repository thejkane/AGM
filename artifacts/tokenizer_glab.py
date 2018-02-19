from __future__ import division
import sys
import re
import glob

print "========== Starting Tokenizer =============="
data = []
header = ('Nodes', 'Cores', 'Binary', 'Threads', 'Scale', 'Sources', 'TMEAN','TSTDEV')
for inputfile in glob.glob("*.o*"):
	lasttuple = ()
	currentscalenodetpl = ()
	with open(inputfile) as input:
	    for line in input:
		#Extract scale and nodes
#                              aprun -n 1 -N 1 -d 32 -b -e GRAPHLAB_SUBNET_ID=10.128.0.0 -e GRAPHLAB_SUBNET_MASK=255.252.0.0 -e GRAPHLAB_THREADS_PER_WORKER=1 ./sssp --rmat 19 --numsources 8 --max_degree_source 1 --rmatversion 1
#aprun -n 4 -N 1 -d 32 -b -e GRAPHLAB_SUBNET_ID=10.128.0.0 -e GRAPHLAB_SUBNET_MASK=255.252.0.0 -e GRAPHLAB_THREADS_PER_WORKER=32 ./sssp --rmat 26 --numsources 4 --max_degree_source 1 --rmatversion 1

		m = re.search(r'aprun -n ([0-9]+) -N ([0-9]+) -d ([0-9]+) -b -e GRAPHLAB_SUBNET_ID=10.128.0.0 -e GRAPHLAB_SUBNET_MASK=255.252.0.0 -e GRAPHLAB_THREADS_PER_WORKER=([0-9]+) ./(\w+) --rmat ([0-9]+) --numsources ([0-9]+) --max_degree_source 1 --rmatversion 1', line)
		if m != None:
		    tasks = int(m.group(1))
		    taskpernode = int(m.group(2))
		    systhreads = int(m.group(3))
		    threads = int(m.group(4))
		    binary = m.group(5)
		    scale = m.group(6)
		    sources = m.group(7)
		    nodes = tasks / taskpernode
		    cores = tasks * threads
		    currentscalenodetpl = (nodes, cores, binary, threads, scale, sources)


		# extract times
		m = re.search(r'MEAN : ([0-9]+.[0-9]+) STDDEV : ([0-9]+.[0-9]+) MIN : ([0-9]+.[0-9]+) Q1 : ([0-9]+.[0-9]+) MEDIAN : ([0-9]+.[0-9]+) Q3 : ([0-9]+.[0-9]+) MAX : ([0-9]+.[0-9]+)', line)
		m = re.search(r'MEAN : ([0-9]+.[0-9]+) STDDEV : ([0-9]+.[0-9]+)', line)
		if m != None:
		    timetple = (m.group(1), m.group(2))
		    lasttuple = currentscalenodetpl + timetple
		    data.append(lasttuple)
		    lasttuple = ()

import csv
outputfile =  "results.csv"
print 'Writing output to ' + outputfile
with open(outputfile, 'wb') as csvfile:
    resultwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    resultwriter.writerow(header)
    for d in data:
        resultwriter.writerow(d)


