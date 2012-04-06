import sys
import math
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

    # read in all the arguments
    for key, val in optlist:
        if key == '-c':
            k = int(val)
        elif key == '-v':
            var_cutoff = float(val)
        elif key == '-i':
            inputname = val

    # check requierd arguments are inputted and valid
    if k < 0 or var_cutoff < 0 or \
            inputname is None:
        usage()
        sys.exit()
    return (k, var_cutoff, inputname)

'''
A wrapper class for a dna strand.
Purpose is to allow strands to be inserted into a min heap when generating
the initial centroids.
'''
class heapStrand:
    def __init__(self, dnastrand, value):
        self.dnastrand = dnastrand
        self.value = value
    def __cmp__(self, obj):
        return cmp(self.value, obj.value)

'''
A class to represent Clusters of dna strands. each cluster contains centroids
and dna strands belonging to the cluster
'''
class Cluster:
    def __init__(self, dnastrands):
        if len(dnastrands) == 0: raise Exception("empty cluster")
        self.dnastrands = dnastrands
        self.dimension = len(dnastrands[0])
        for d in dnastrands:
            if len(d) != self.dimension:
		raise Exception("strands not all of same length")
        self.centroid = self.calculateCentroid()
    def __repr__(self):
        return str(self.dnastrands)
    # when points are added to clusters, centroids need to be updated
    def update(self, dnastrands):
        old_centroid = self.centroid
        self.dnastrands = dnastrands
        self.centroid = self.calculateCentroid()
        return getDistance(old_centroid, self.centroid)
    # calculates the centroid of the cluster based on all data points in it
    def calculateCentroid(self):
        centroid_strand = []
	# iterates through the bases of each dna strand
        for i in range (len(self.dnastrands[0])):
            aCount, tCount, cCount, gCount = 0, 0, 0, 0
            # gather counts of bases at that index
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
            # find which base has max count
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

'''
Runs the k-means clustering algorithm given k and a set of dna strands..
i.e. partitions the dna strands into k clusters based on similarity
to the centroids of the clusters.
'''
def kmeans(dnastrands, k, var_cutoff, dnaLength):

    #Step 1: generate an initial array of centroid dna strands

    strandsSet = Set(dnastrands)
    initial = []
    for i in range(k):
	# randomly select a dna strand to be a centroid
        centroid = random.choice(list(strandsSet))
        initial.append(centroid)
        # build heap based on distance from chosen centroid
        distHeap = []
        for d in strandsSet:
            hstrand = heapStrand(d, getDistance(centroid, d))
            heapq.heappush(distHeap, hstrand)
	# remove (#points/k) points from the data set so that strands closest
	# to the already chosen centroids are not also chosen as centroids
        for j in range(len(dnastrands)/k):
            toRemove = heapq.heappop(distHeap).dnastrand
            strandsSet.remove(toRemove)
    centroids = [Cluster([d]) for d in initial]

    # Step 2: run k means on initial centroids until cutoff_var is reached

    while True:
        clusters = [[] for c in centroids]
        for d in dnastrands:
            index = 0
	    # calculate the distance from that dna to all strands
            smallest_distance = float('inf')
            for i in range(len(centroids)):
                distance = getDistance(d, centroids[i].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i
	    # the strand belongs to the centroid of least distance to it
            clusters[index].append(d)
        max_centroidvar = float('-inf')
	# calculates the centroid and updates variance
        for i in range(len(centroids)):
            var = centroids[i].update(clusters[i])
            max_centroidvar = max(max_centroidvar, var)
	# algo is finished when variance threshold is satisfied in all clusters
        if max_centroidvar < var_cutoff:
            break
    return centroids

'''
Returns the Hamming Distance between two given strands
'''
def getDistance(s1, s2):
    if len(s1) != len(s2): raise Exception("dnastrands are not the same length")
    diff = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            diff += 1
    return diff

def main():
    k, var_cutoff, inputname = handleArgs(sys.argv)
    # read dna strands from input file and add them to the dnastrands array
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
    # run the k-means algorithm on the dna strands
    clusters = kmeans(dnastrands, k, var_cutoff, dnaLength)

    for i,c in enumerate(clusters):
        count = 0
        for d in c.dnastrands:
            print " Cluster: ",i+1,"\t DNA Strand:", d
            count += 1
        print 'cluster ' + str(i+1) +': ' + str(count) + ' Strands'
if __name__ == "__main__":
    cProfile.run("main()")
