import sys
import math
import numpy
import random
import getopt

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
        centroid_strand = []
        for i in len(dnastrands[0]):
            aCount, tCount, cCount, gCount = 0, 0, 0, 0
            # gather counts
            for d in dnastrands:
                if d[i] == 'a':
                    aCount += 1
                elif d[i] == 't':
                    tCount += 1
                elif d[i] == 'c':
                    cCount += 1
                elif d[i] == 'g':
                    gCount += 1
            maxCount = max([aCount, tCount, cCount, gCount])
            # find which had max count
            mostCommon = []
            if aCount == maxCount:
                mostCommon.append('a')
            elif tCount == maxCount:
                mostCommon.append('t')
            elif cCount == maxCount:
                mostCommon.append('c')
            elif gCount == maxCount:
                mostCommon.append('g')
            if len(mostCommon) > 1:
                # tie exists, so randomly choose one
                index = random.randint(0,len(mostCommon)-1)
                centroid_strand.append(mostCommon[index])
            else:
                centroid_strand.append(mostCommon[0])
        return tuple(centroid_strand)

def kmeans(dnastrands, k, var_cutoff):
    """choose k initial centroids randomly, but check
    similarity and regenerate if too similar"""
    # DIXIE: write this
    initial = []
    min_coord = map(min,zip(*points))
    max_coord = map(max,zip(*points))
    dim = len(max_coord)
    for n in range(1,k+1):
        c = []
        for d in range(0,dim):
            p = (max_coord[d]-min_coord[d])/k*n
            c.append(p)
        initial.append(tuple(c))
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

def getDistance(s1, s2):
    # find hamming distance between s1 and s2
    if len(s1) != len(s2): raise Exception("dnastrands are not the same length")
    diff = 0
    for i in len(s1):
        if s1[i] != s2[i]:
            diff += 1
    return diff

def main():
    # TODO: rewrite this to read dnastrand input
    k, var_cutoff, inputname = handleArgs(sys.argv)
    # read points from input file
    inputfile = open(inputname, "rb")
    inputstream = csv.reader(inputfile)
    rownum = 0
    points = []
    for row in inputstream:
        colnum = 0
        for col in row:
            if colnum == 0:
                x = float(col)
            if colnum == 1:
                y = float(col)
            colnum += 1
        point = tuple([x,y])
        points.append(point)
    inputfile.close()
    # end
    clusters = kmeans(dnastrands, k, var_cutoff)

    for i,c in enumerate(clusters):
        for d in c.dnastrands:
            print " Cluster: ",i,"\t DNA Strand:", d

if __name__ == "__main__":
    main()