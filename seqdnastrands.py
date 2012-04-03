import sys
import math
import numpy
import random
import getopt
import csv

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
        # TONY: write this
        reduce_coord = lambda i:reduce(lambda x,p : x + p[i],self.points,0.0)
        centroid_coords = [reduce_coord(i)/len(self.points) for i in range(self.dimension)]
        return tuple(centroid_coords)
        # TONY: end

#creates a centroid with a given length, by randomly selecting bases
def createCentroid(dnaLength):
   dna = []
   for j in range(dnaLength):
       dna.append(random.choice(bases))
   return dna

#returns the similarity between two given centroids: 'centroid' & 'pair'
def similarity(centroid, pair, dnaLength):
   sameCount = 0
   for i in range(dnaLength):
      if centroid[i] == pair[i]:
         sameCount += 1
   return sameCount

#returns the most similarities between a given centroid and
#all other existing centroids
def maxSimilarity(centroid, initial, dnaLength):
   maxSim = 0
   for pair in initial:
      sim = similarity(centroid, pair, dnaLength)
      if sim > maxSim:
         maxSim = sim
   return maxSim

def kmeans(dnastrands, k, var_cutoff, dnaLength):
    """choose k initial centroids randomly, but check
    similarity and regenerate if too similar"""
    # DIXIE: write this
    initial = []
    maxAllowedSim = int(dnaLength * .75)
    for i in range(k):
        centroid = createCentroid(dnaLength)
        tryCount = 0
        #is it different enough from other centroids in initial?
        similarity = maxSimilarity(centroid, initial, dnaLength)
        minSimSoFar = similarity
        bestCentroid = centroid
        while (similarity > maxAllowedSim):
            #if we've tried too many times to satisfy maxAllowedSim,
            #just take the most different centroid we've generated
            if (tryCount == 20):
                break
            if (similarity < minSimSoFar):
                minSimSoFar = similarity
                bestCentroid = centroid
            centroid = createCentroid(dnaLength)
            similarity = maxSimilarity(centroid, initial, dnaLength)
            tryCount += 1
        initial.append(bestCentroid)

    # DIXIE: end
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

def getDistance(x, y):
    # TONY: edit this, make it hamming distance
    if len(x) != len(y): raise Exception("coordinates not in same dimension")
    ret = reduce(lambda a,b: a + pow((x[b]-y[b]), 2), range(len(x)), 0.0)
    return math.sqrt(ret)
    # TONY: end

def main():
    # TODO: rewrite this to read dnastrand input
    k, var_cutoff, inputname = handleArgs(sys.argv)
    # read points from input file
    inputfile = open(inputname, "rb")
    inputstream = csv.reader(inputfile)
    strands = []
    for strand in inputstream:
        strands.append(strand[0])
    inputfile.close()
    dnaLength = len(strands[0])
    # end
    clusters = kmeans(dnastrands, k, var_cutoff, dnaLength)

    for i,c in enumerate(clusters):
        for d in c.dnastrands:
            print " Cluster: ",i,"\t DNA Strand:", d

if __name__ == "__main__":
    main()
