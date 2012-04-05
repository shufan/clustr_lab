import sys
import math
import numpy
import random
import getopt
import csv
import heapq
import cProfile
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

# #creates a centroid with a given length, by randomly selecting bases
# def createCentroid(dnaLength):
#    dna = []
#    for j in range(dnaLength):
#        dna.append(random.choice(bases))
#    return dna

# #returns the similarity between two given centroids: 'centroid' & 'pair'
# def similarity(centroid, pair, dnaLength):
#    sameCount = 0
#    for i in range(dnaLength):
#       if centroid[i] == pair[i]:
#          sameCount += 1
#    return sameCount

# #returns the most similarities between a given centroid and
# #all other existing centroids
# def maxSimilarity(centroid, initial, dnaLength):
#    maxSim = 0
#    for pair in initial:
#       sim = similarity(centroid, pair, dnaLength)
#       if sim > maxSim:
#          maxSim = sim
#    return maxSim

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
    while True:
        clusters = [[] for c in centroids]
        for d in dnastrands:
            index = 0
            smallest_distance = float('inf')
            for i in range(len(centroids)):
                distance = getDistance(d, centroids[i].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i
            clusters[index].append(d)
        max_centroidvar = float('-inf')
        for i in range(len(centroids)):
            var = centroids[i].update(clusters[i])
            max_centroidvar = max(max_centroidvar, var)
        if max_centroidvar < var_cutoff:
            break
    return centroids

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
