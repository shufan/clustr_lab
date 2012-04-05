import sys
import math
import numpy
import random
import getopt
import csv
import heapq
import cProfile
from mpi4py import MPI
from sets import Set

bases = ['A', 'C', 'G', 'T']

def usage():
    print '$> python seqdatapoints.py <required args> [optional args]\n' + \
        '\t-c <#>\t\tNumber of clusters to find\n' + \
        '\t-v <#>\t\tCentroid variance cutoff\n' + \
        '\t-i <file>\tFilename for the output of the raw data\n'

def handleArgs(args):
    #set up return values
    k = -1
    var_cutoff = -1
    inputname = None

    try:
        optlist, args = getopt.getopt(args[1:], 'c:v:i:')
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(2)

    for key, val in optlist:
        if key == '-c':
            k = int(val)
        elif key == '-v':
            var_cutoff = float(val)
        elif key == '-i':
            inputname = val

    if k < 0 or var_cutoff < 0 or \
            inputname is None:
        usage()
        sys.exit()
    return (k, var_cutoff, inputname)

class heapStrand:
    def __init__(self, dnastrand, value):
        self.dnastrand = dnastrand
        self.value = value
    def __cmp__(self, obj):
        return cmp(self.value, obj.value)

class Cluster:
    def __init__(self, dnastrands):
        if len(dnastrands) == 0: raise Exception("empty cluster")
        self.dnastrands = dnastrands
        self.dimension = len(dnastrands[0])
        for d in dnastrands:
            if len(d) != self.dimension: raise Exception("strands not all of same length")
        self.centroid = self.calculateCentroid()
    def __repr__(self):
        return str(self.dnastrands)
    def update(self, dnastrands):
        old_centroid = self.centroid
        self.dnastrands = dnastrands
        self.centroid = self.calculateCentroid()
        return getDistance(old_centroid, self.centroid)
    def calculateCentroid(self):
        centroid_strand = []
        for i in range (len(self.dnastrands[0])):
            aCount, tCount, cCount, gCount = 0, 0, 0, 0
            # gather counts
            for d in self.dnastrands:
                if d[i] == 'A':
                    aCount += 1
                elif d[i] == 'T':
                    tCount += 1
                elif d[i] == 'C':
                    cCount += 1
                elif d[i] == 'G':
                    gCount += 1
            maxCount = max([aCount, tCount, cCount, gCount])
            # find which had max count
            mostCommon = []
            if aCount == maxCount:
                mostCommon.append('A')
            elif tCount == maxCount:
                mostCommon.append('T')
            elif cCount == maxCount:
                mostCommon.append('C')
            elif gCount == maxCount:
                mostCommon.append('G')
            if len(mostCommon) > 1:
                # tie exists, so randomly choose one
                index = random.randint(len(mostCommon))
                centroid_strand.append(mostCommon[index])
            else:
                centroid_strand.append(mostCommon[0])
        return tuple(centroid_strand)

def partition(lst, n, needIndices=False):
    ''' partitioning code from http://stackoverflow.com/questions/2659900/
    python-slicing-a-list-into-n-nearly-equal-length-partitions'''
    division = len(lst) / float(n)
    chunkList = [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n)]
    if needIndices:
        chunkIndexList = []
        for i in xrange(0,n):
            chunkIndexList.append(int(round(i * division)))
        return (chunkList, chunkIndexList)
    return chunkList

def kmeans(dnastrands, k, var_cutoff, dnaLength):
    strandsSet = Set(dnastrands)
    initial = []
    for i in range(k):
        centroid = random.choice(list(strandsSet))
        initial.append(centroid)
        #build heap based on distance from chosen centroid
        distHeap = []
        for d in strandsSet:
            hstrand = heapStrand(d, getDistance(centroid, d))
            heapq.heappush(distHeap, hstrand)
        for j in range(len(dnastrands)/k):
            toRemove = heapq.heappop(distHeap).dnastrand
            strandsSet.remove(toRemove)
    # run k means on initial centroids until cutoff_var reached
    centroids = [Cluster([d]) for d in initial]
    # mpi communicator
    comm = MPI.COMM_WORLD

    mpirank = comm.Get_rank()
    mpisize = comm.Get_size()

    while True:
        clusters = [c for c in centroids]
        clusters = comm.bcast(clusters, root=0)
        scatteredStrands = comm.scatter(partition(dnastrands, mpisize), root=0)
        chunkedStrandClstrMap = {}
        for d in scatteredStrands:
            index = 0
            smallest_distance = float('inf')
            for i in range(len(clusters)):
                distance = getDistance(d, clusters[i].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i
                    # if index == 0:
                    #     print 'HAS 0\n'
                    # elif index == 1:
                    #     print 'HAS 1\n'
                    # elif index == 2:
                    #     print 'HAS 2\n'
                    # elif index == 3:
                    #     print 'HAS 3\n'
                    # elif index == 4:
                    #     print 'HAS 4\n'
            chunkedStrandClstrMap.setdefault(index, []).append(d)
            # print 'this should never be 0: count for cluster ' + str(index) + ' --->' + str(len(chunkedStrandClstrMap[index]))
        # for i in range(len(centroids)):
            # print 'number in cluster' + str(i) + 'for this chunk' + str(len(chunkedStrandClstrMap[i]))
        chunkedStrandClstrMap = comm.gather(chunkedStrandClstrMap, root=0)

        strandClstrMap = {}
        if mpirank == 0:
            for m in chunkedStrandClstrMap:
                for key in m.keys():
                    if strandClstrMap.setdefault(key, Set(m.get(key, []))) != Set(m.get(key, [])):
                        strandClstrMap[key].update(Set(m.get(key, [])))
            # print 'cluster mappings lengths...' + str(len(strandClstrMap[4]))
        strandClstrMap = comm.bcast(strandClstrMap, root=0)

        (chunkedClusters, procIndices) = partition(clusters, mpisize, True)
        scatteredClusters = comm.scatter(chunkedClusters, root=0)
        procIndices = comm.bcast(procIndices, root=0)
        for i in range(len(scatteredClusters)):
            scatteredClusters[i] = []
            [scatteredClusters[i].append(d) for d in strandClstrMap.get(procIndices[mpirank]+i, [])]

        gatheredClusters = comm.gather(scatteredClusters, root=0)
        if mpirank == 0:
            gatheredClusters = reduce(lambda x, y: x+y, gatheredClusters)

        done = False

        if mpirank == 0:
            max_centroidvar = float('-inf')
            for i in range(len(centroids)):
                # print 'updating with...' + str(gatheredClusters[i]) + 'on iteration' + str(i)
                var = centroids[i].update(gatheredClusters[i])
                max_centroidvar = max(max_centroidvar, var)
            if max_centroidvar < var_cutoff:
                done = True

        done = comm.bcast(done, root=0)
        if done:
            break
    if mpirank == 0:
        return centroids
    else:
        return []

def getDistance(s1, s2):
    # find hamming distance between s1 and s2
    if len(s1) != len(s2): raise Exception("dnastrands are not the same length")
    diff = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    return diff

def main():
    k, var_cutoff, inputname = handleArgs(sys.argv)
    # read points from input file
    inputfile = open(inputname, "rb")
    inputstream = csv.reader(inputfile)
    dnastrands = []
    for d in inputstream:
        dlist = []
        for i in range(len(d[0])):
            dlist.append(d[0][i])
        dnastrands.append(tuple(dlist))
    inputfile.close()
    dnaLength = len(dnastrands[0])
    clusters = kmeans(dnastrands, k, var_cutoff, dnaLength)

    for i,c in enumerate(clusters):
        count = 0
        for d in c.dnastrands:
            print " Cluster: ",i,"\t DNA Strand:", d
            count += 1
        print 'cluster ' + str(i) +': ' + str(count)
if __name__ == "__main__":
    cProfile.run("main()")
