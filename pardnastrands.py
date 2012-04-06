import sys
import math
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
A class to represent Clusters of dna strands. each cluster contains
dna strands belonging to the cluster
'''
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
        # when strands are added to clusters, centroids need to be updated
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

def partition(lst, n, needIndices=False):
    ''' partitioning code from http://stackoverflow.com/questions/2659900/
    python-slicing-a-list-into-n-nearly-equal-length-partitions
    description of function in the url above'''
    division = len(lst) / float(n)
    chunkList = [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n)]
    if needIndices:
        chunkIndexList = []
        for i in xrange(0,n):
            chunkIndexList.append(int(round(i * division)))
        return (chunkList, chunkIndexList)
    return chunkList

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

    # mpi communicator
    comm = MPI.COMM_WORLD

    mpirank = comm.Get_rank()
    mpisize = comm.Get_size()

    while True:
        clusters = [c for c in centroids]
        # broadcast clusters and scatter the data among processors
        clusters = comm.bcast(clusters, root=0)
        scatteredStrands = comm.scatter(partition(dnastrands, mpisize), root=0)
        # initialize empty map of clusters to the strands nearest to them
        chunkedStrandClstrMap = {}
        # build a map on each machine
        for d in scatteredStrands:
            index = 0
            # calculate the distance from that dna to all strands
            smallest_distance = float('inf')
            for i in range(len(clusters)):
                distance = getDistance(d, clusters[i].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i
            # the strand belongs to the centroid of least distance to it
            chunkedStrandClstrMap.setdefault(index, []).append(d)
        # gather all maps from each machine into a list of the maps
        chunkedStrandClstrMap = comm.gather(chunkedStrandClstrMap, root=0)
        # combine the information on the list of maps into one unified map
        strandClstrMap = {}
        if mpirank == 0:
            for m in chunkedStrandClstrMap:
                for key in m.keys():
                    if strandClstrMap.setdefault(key, Set(m.get(key, []))) != Set(m.get(key, [])):
                        strandClstrMap[key].update(Set(m.get(key, [])))
        # broadcast the map of clusters and their strands to all the machines
        strandClstrMap = comm.bcast(strandClstrMap, root=0)
        # partition the clusters and scatter them among the machines
        (chunkedClusters, procIndices) = partition(clusters, mpisize, True)
        scatteredClusters = comm.scatter(chunkedClusters, root=0)
        # broadcast the original indices of the cluster at the beginning of
        # partition of the clusters in each machine
        procIndices = comm.bcast(procIndices, root=0)
        # add the proper strands for each cluster to the cluster
        for i in range(len(scatteredClusters)):
            scatteredClusters[i] = []
            [scatteredClusters[i].append(d) for d in strandClstrMap.get(procIndices[mpirank]+i, [])]

        # gather the updated clusters into a list of lists of clusters
        gatheredClusters = comm.gather(scatteredClusters, root=0)
        if mpirank == 0
            # reduce it to a list of clusters
            gatheredClusters = reduce(lambda x, y: x+y, gatheredClusters)

        done = False

        # master updates clusters to the gathered clusters, and gets centroids
        if mpirank == 0:
            max_centroidvar = float('-inf')
            # calculates the centroid and updates variance
            for i in range(len(centroids)):
                var = centroids[i].update(gatheredClusters[i])
                max_centroidvar = max(max_centroidvar, var)
            # algo is finished when variance threshold is satisfied in all clusters
            if max_centroidvar < var_cutoff:
                done = True
        # broadcast to all machines that the algorithm is done
        done = comm.bcast(done, root=0)
        if done:
            break
    if mpirank == 0:
        # master returns the calculated centroids
        return centroids
    else:
        return []

'''
Returns the Hamming Distance between two given strands
'''
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
