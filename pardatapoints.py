# Tony Zhang (tszhang) & Dixie Kee (dkee)

import sys
import math
import csv
import random
import getopt
import heapq
import cProfile
from mpi4py import MPI
from sets import Set

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

    # read in all the arguments
    for key, val in optlist:
        if key == '-c':
            k = int(val)
        elif key == '-v':
            var_cutoff = float(val)
        elif key == '-i':
            inputname = val

    # check required arguments are inputted and valid
    if k < 0 or var_cutoff < 0 or \
            inputname is None:
        usage()
        sys.exit()
    return (k, var_cutoff, inputname)

'''
A wrapper class for a point.
Purpose is to allow points to be inserted into a min heap when generating
the intial centroids.
'''
class heapPoint:
    def __init__(self, point, value):
        self.point = point
        self.value = value
    def __cmp__(self, obj):
        return cmp(self.value, obj.value)

'''
A class to represent Clusters of data points. each cluster contains
data points belonging to the cluster
'''
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
    # when points are added to clusters, centroids need to be updated
    def update(self, points):
        old_centroid = self.centroid
        self.points = points
        self.centroid = self.calculateCentroid()
        return getDistance(old_centroid, self.centroid)
    # calculates the centroid of the cluster based on all data points in it
    def calculateCentroid(self):
        reduce_coord = lambda i:reduce(lambda x,p : x + p[i],self.points,0.0)
        centroid_coords = [reduce_coord(i)/len(self.points) \
						for i in range(self.dimension)]
        return tuple(centroid_coords)

def partition(lst, n, needIndices=False):
    ''' partitioning code from http://stackoverflow.com/questions/2659900/
    python-slicing-a-list-into-n-nearly-equal-length-partitions
    description of function in the url above'''
    division = len(lst) / float(n)
    chunkList = [lst[int(round(division * i)): \
			int(round(division * (i + 1)))] for i in xrange(n)]
    if needIndices:
        chunkIndexList = []
        for i in xrange(0,n):
            chunkIndexList.append(int(round(i * division)))
        return (chunkList, chunkIndexList)
    return chunkList

'''
Runs the k-means clustering algorithm given k and a set of datapoints.
i.e. partitions the data points into k clusters based on similarity
to the centroids of the clusters.
'''
def kmeans(points, k, var_cutoff):
    # Step 1: generate an initial array of centroid data points
    pointsSet = Set(points)
    initial = []
    for i in range(k):
        # randomly select a data point to be a centroid
        centroid = random.choice(list(pointsSet))
        initial.append(centroid)
        # build a min heap based on distance from chosen centroid
        distHeap = []
        for p in pointsSet:
            hpoint = heapPoint(p,getDistance(centroid, p))
            heapq.heappush(distHeap, hpoint)
        # remove (#points/k) points from the data set so that points closest
        # to the already chosen centroid are not also chosen as centroids
        for j in range(len(points)/k):
            toRemove = heapq.heappop(distHeap).point
            pointsSet.remove(toRemove)
    centroids = [Cluster([p]) for p in initial]

    # Step 2: run k means on initial centroids until cutoff_var is reached

    # mpi communicator
    comm = MPI.COMM_WORLD

    mpirank = comm.Get_rank()
    mpisize = comm.Get_size()

    while True:
        clusters = [c for c in centroids]
        # broadcast clusters and scatter the data among processors
        clusters = comm.bcast(clusters, root=0)
        scatteredPoints = comm.scatter(partition(points, mpisize), root=0)
        # initialize empty map of clusters to the points nearest to them
        chunkedPtClstrMap = {}
        # build a map on each machine
        for p in scatteredPoints:
            index = 0
            # calculate the distance from that point to all points
            smallest_distance = float('inf')
            for i in range(len(clusters)):
                distance = getDistance(p, clusters[i].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i
            chunkedPtClstrMap.setdefault(index, []).append(p)
        # gather all maps from each machine into a list of the maps
        chunkedPtClstrMap = comm.gather(chunkedPtClstrMap, root=0)
        # combine the information on the list of maps into one unified map
        ptClstrMap = {}
        if mpirank == 0:
            for m in chunkedPtClstrMap:
                for key in m.keys():
                    if ptClstrMap.setdefault(key, Set(m.get(key, []))) != \
							Set(m.get(key, [])):
                        ptClstrMap[key].update(Set(m.get(key, [])))
        # broadcast the map of clusters and their points to all the machines
        ptClstrMap = comm.bcast(ptClstrMap, root=0)
        # partition the clusters and scatter them among the machines
        (chunkedClusters, procIndices) = partition(clusters, mpisize, True)
        scatteredClusters = comm.scatter(chunkedClusters, root=0)
        # broadcast the original indices of the cluster at the beginning of
        # partition of the clusters in each machine
        procIndices = comm.bcast(procIndices, root=0)
        # add the proper points for each cluster to the cluster
        for i in range(len(scatteredClusters)):
            scatteredClusters[i] = []
            [scatteredClusters[i].append(p) \
			for p in ptClstrMap.get(procIndices[mpirank]+i, [])]

        # gather the updated clusters into a list of lists of clusters
        gatheredClusters = comm.gather(scatteredClusters, root=0)
        if mpirank == 0:
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
            # algo is finished when variance threshold is satisfied
	    # in all clusters
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
Returns the Euclidean Distance of two points
'''
def getDistance(x, y):
    if len(x) != len(y): raise Exception("coordinates not in same dimension")
    ret = reduce(lambda a,b: a + pow((x[b]-y[b]), 2), range(len(x)), 0.0)
    return math.sqrt(ret)

def main():
    k, var_cutoff, inputname = handleArgs(sys.argv)
    # read points from input file and add them to the points array
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
    # run the k-means algorithm on the points
    clusters = kmeans(points, k, var_cutoff)

    for i,c in enumerate(clusters):
        count = 0
        for p in c.points:
            print " Cluster: ",i+1,"\t Point:", p
            count += 1
        print 'Cluster ' + str(i+1) +': ' + str(count) + ' Points'
if __name__ == "__main__":
    cProfile.run("main()")






