import sys
import math
import numpy
import csv
import random
import getopt

def usage():
    print '$> python seqdatapoints.py <required args> [optional args]\n' + \
        '\t-c <#>\t\tNumber of clusters to find\n' + \
        '\t-v <#>\t\tCentroid variance cutoff\n' + \
        '\t-i <file>\tFilename for the input data\n'

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
    def __init__(self, points):
        if len(points) == 0: raise Exception("empty cluster")
        self.points = points
        self.dimension = len(points[0])
        for p in points:
            if len(p) != self.dimension: raise Exception("illegal coordinate")
        self.centroid = self.calculateCentroid()
    def __repr__(self):
        return str(self.points)
    def update(self, points):
        old_centroid = self.centroid
        self.points = points
        self.centroid = self.calculateCentroid()
        return getDistance(old_centroid, self.centroid)
    def calculateCentroid(self):
        reduce_coord = lambda i:reduce(lambda x,p : x + p[i],self.points,0.0)
        centroid_coords = [reduce_coord(i)/len(self.points) for i in range(self.dimension)]
        return tuple(centroid_coords)

def kmeans(points, k, var_cutoff):
    """choose k initial centroids by taking max
    and min of each dimensino and take k step"""
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
    centroids = [Cluster([p]) for p in initial]
    while True:
        clusters = [[] for c in centroids]
        for p in points:
            index = 0
            smallest_distance = float('inf')
            for i in range(len(centroids)):
                distance = getDistance(p, centroids[i].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i
            clusters[index].append(p)
        max_centroidvar = float('-inf')
        for i in range(len(centroids)):
            var = centroids[i].update(clusters[i])
            max_centroidvar = max(max_centroidvar, var)
        if max_centroidvar < var_cutoff:
            break
    return centroids

def getDistance(x, y):
    if len(x) != len(y): raise Exception("coordinates not in same dimension")
    ret = reduce(lambda a,b: a + pow((x[b]-y[b]), 2), range(len(x)), 0.0)
    return math.sqrt(ret)

def main():
    k, var_cutoff, inputname = handleArgs(sys.argv)
    #read points from input file
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
    clusters = kmeans(points, k, var_cutoff)

    for i,c in enumerate(clusters):
        for p in c.points:
            print " Cluster: ",i,"\t Point:", p

if __name__ == "__main__":
    main()






