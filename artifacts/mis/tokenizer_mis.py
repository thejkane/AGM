from __future__ import division
import sys
import re
import glob

print "========== Starting Tokenizer =============="
data = []
header = ('Nodes', 'Binary', 'Threads', 'Scale', 'Degree', 'Sources', 'Coalescing', 'Flush', 'Algorithm', 'TMEAN','TSTDEV','TMIN','TQ1','TMEDIAN', 'TQ3','TMAX','TEPMIN','TEPQ1','TEPMEDIAN','TEPQ3','TEPMAX','TEPMEAN','TEPSTDEV')
for inputfile in glob.glob("*.o*"):
	lasttuple = ()
	currentscalenodetpl = ()
	with open(inputfile) as input:
	    for line in input:
		#Extract scale and nodes
		m = re.search(r'\+ aprun -b -n ([0-9]+) -N ([0-9]+) -d ([0-9]+) -r 1 ./(\w+) --poll-task 1 --threads ([0-9]+) --scale ([0-9]+) --degree ([0-9]+) --num-sources ([0-9]+)', line)
		if m != None:
		    tasks = int(m.group(1))
		    taskpernode = int(m.group(2))
		    nodes = tasks / taskpernode
		    currentscalenodetpl = (nodes, m.group(4), m.group(5), m.group(6), m.group(7), m.group(8))

		m = re.search(r'\+ aprun -b -n ([0-9]+) -N ([0-9]+) -d ([0-9]+) ./(\w+) --poll-task 1 --threads ([0-9]+) --scale ([0-9]+) --degree ([0-9]+) --num-sources ([0-9]+)', line)
		if m != None:
		    tasks = int(m.group(1))
		    taskpernode = int(m.group(2))
		    nodes = tasks / taskpernode
		    currentscalenodetpl = (nodes, m.group(4), m.group(5), m.group(6), m.group(7), m.group(8))


		m = re.search(r'\[(\w+)\-(\w+)\] Data Structure :(\w+) Threads: ([0-9]+) Coalescing: ([0-9]+) Poll: ([0-9]+) Routing: ([0-9]+) Depth: ([0-9]+) Flush: ([0-9]+)per_thread_reductions - ([0-9]+) no_reductions : ([0-9]+)', line)
		if m != None:
		    algo=m.group(1)+'-'+m.group(2)
		    coals=m.group(5)
		    flsh=m.group(9)
		    lasttuple = currentscalenodetpl + (coals,flsh,algo)

		m = re.search(r'\[(\w+)\-(\w+)\] Luby Algorithm :(\w+) Data Structure :(\w+) Threads: ([0-9]+) Coalescing: ([0-9]+) Poll: ([0-9]+) Routing: ([0-9]+) Depth: ([0-9]+) Flush: ([0-9]+)per_thread_reductions - ([0-9]+) no_reductions : ([0-9]+)', line)
		if m != None:
		    algo=m.group(1)+'-'+m.group(2)+'-'+m.group(3)
		    coals=m.group(6)
		    flsh=m.group(10)
		    lasttuple = currentscalenodetpl + (coals,flsh,algo)

		# extract times
		m = re.search(r'MEAN : ([0-9]+.[0-9]+) STDDEV : ([0-9]+.[0-9]+) MIN : ([0-9]+.[0-9]+) Q1 : ([0-9]+.[0-9]+) MEDIAN : ([0-9]+.[0-9]+) Q3 : ([0-9]+.[0-9]+) MAX : ([0-9]+.[0-9]+)', line)
		if m != None:
		    timetple = (m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6), m.group(7))
		    lasttuple = lasttuple + timetple
		    data.append(lasttuple)
		    lasttuple = ()

		# extract TEPS
#		m = re.search(r'Min\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
#		if m != None:
#		    minteptpl = (m.group(1),)
#		    lasttuple = lasttuple + minteptpl
#		m = re.search(r'First Quartile\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
#		if m != None:
#		    q1teptpl = (m.group(1),)
#		    lasttuple = lasttuple + q1teptpl
#		m = re.search(r'Median\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
#		if m != None:
#		    medteptpl = (m.group(1),)
#		    lasttuple = lasttuple + medteptpl
#		m = re.search(r'Third Quartile\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
#		if m != None:
#		    q3teptpl = (m.group(1),)
#		    lasttuple = lasttuple + q3teptpl
#		    m = re.search(r'Max\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
#		if m != None:
#		    maxteptpl = (m.group(1),)
#		    lasttuple = lasttuple + maxteptpl
#		m = re.search(r'Harmonic Mean\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
#		if m != None:
#		    meanteptpl = (m.group(1),)
#		    lasttuple = lasttuple + meanteptpl
#		m = re.search(r'Harmonic Stddev\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
#		if m != None:
#		    stdteptpl = (m.group(1),)
#		    lasttuple = lasttuple + stdteptpl
#		    data.append(lasttuple)
#		    lasttuple = ()
		    
import csv
outputfile =  "results.csv"
print 'Writing output to ' + outputfile
with open(outputfile, 'wb') as csvfile:
    resultwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    resultwriter.writerow(header)
    for d in data:
        resultwriter.writerow(d)


