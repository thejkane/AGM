from __future__ import division
import sys
import re
inputfile = sys.argv[1]

if inputfile == None:
    print "[Error] Enter valid file name !"


print "========== Starting Tokenizer =============="
data = []
header = ('Scale', 'Nodes', 'Algorithm', 'DS', 'Threads', 'Coalescing', 'Poll', 'Routing','Depth','Delta','Sources','Useful','Invalid', 'DJRatio','Useless','Rejected','TMEAN','TSTDEV','TMIN','TQ1','TMEDIAN', 'TQ3','TMAX','TEPMIN','TEPQ1','TEPMEDIAN','TEPQ3','TEPMAX','TEPMEAN','TEPSTDEV')
lasttuple = ()
currentscalenodetpl = ()
with open(inputfile) as input:
    for line in input:
        #Extract scale and nodes
        m = re.search(r'Running Scale : ([0-9]+), Nodes : ([0-9]+)', line)
        if m != None:
            currentscalenodetpl = (m.group(1), m.group(2))

        m = re.search(r'\[(\w+)\] Data Structure :(\w+) Threads: ([0-9]+) Coalescing: ([0-9]+) Poll: ([0-9]+) Routing: ([0-9]+) Depth: ([0-9]+) Delta: ([0-9]+)', line)
        if m != None:
            lasttuple = currentscalenodetpl + (m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6), m.group(7), m.group(8))
        # extract number of sources
        m = re.search('Total [A-Za-z0-9_\- \(\)]+ time for ([0-9]+) sources', line)
        if m != None:
            sourcetple = (m.group(1),)
            lasttuple = lasttuple + sourcetple

        # extract work stats
        m = re.search(r'Useful work: ([0-9]+) \(per source\), invalidated work: ([0-9]+) \(per source\), useless work: ([0-9]+) \(per source\), rejected work: ([0-9]+) \(per source\).', line)
        # Useful, Invalidated, Useless, Rejected
        if m != None:
            useful = int(m.group(1))
            invalidated = int(m.group(2))
            djratio = useful / (useful - invalidated)
            wrktuple = (useful, invalidated, djratio, m.group(3), m.group(4))
            lasttuple = lasttuple + wrktuple

        # extract times
        m = re.search(r'MEAN : ([0-9]+.[0-9]+) STDDEV : ([0-9]+.[0-9]+) MIN : ([0-9]+.[0-9]+) Q1 : ([0-9]+.[0-9]+) MEDIAN : ([0-9]+.[0-9]+) Q3 : ([0-9]+.[0-9]+) MAX : ([0-9]+.[0-9]+)', line)
        if m != None:
            timetple = (m.group(1), m.group(2), m.group(3), m.group(4), m.group(5), m.group(6), m.group(7))
            lasttuple = lasttuple + timetple

        # extract TEPS
        m = re.search(r'Min\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
        if m != None:
            minteptpl = (m.group(1),)
            lasttuple = lasttuple + minteptpl
        m = re.search(r'First Quartile\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
        if m != None:
            q1teptpl = (m.group(1),)
            lasttuple = lasttuple + q1teptpl
        m = re.search(r'Median\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
        if m != None:
            medteptpl = (m.group(1),)
            lasttuple = lasttuple + medteptpl
        m = re.search(r'Third Quartile\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
        if m != None:
            q3teptpl = (m.group(1),)
            lasttuple = lasttuple + q3teptpl
        m = re.search(r'Max\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
        if m != None:
            maxteptpl = (m.group(1),)
            lasttuple = lasttuple + maxteptpl
        m = re.search(r'Harmonic Mean\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
        if m != None:
            meanteptpl = (m.group(1),)
            lasttuple = lasttuple + meanteptpl
        m = re.search(r'Harmonic Stddev\(TEPS\): ([0-9]+.[0-9]+e\+[0-9]+)', line)
        if m != None:
            stdteptpl = (m.group(1),)
            lasttuple = lasttuple + stdteptpl
            data.append(lasttuple)
            lasttuple = ()

import csv
outputfile = inputfile + ".csv"
print 'Writing output to ' + outputfile
with open(outputfile, 'wb') as csvfile:
    resultwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    resultwriter.writerow(header)
    for d in data:
        resultwriter.writerow(d)


