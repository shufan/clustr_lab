import sys
import math
import csv
import random
import getopt
import heapq
import cProfile
from sets import Set

def usage():
    print '$> python seqdatapoints.py <required args> [optional args]\n' + \
        '\t-c <#>\t\tNumber of clusters to find\n' + \
        '\t-v <#>\t\tCentroid variance cutoff\n' + \
        '\t-i <file>\tFilename for the input data\n'

def handleArgs(args):
    # set up return values
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

    # check required aruguments are inputted and valid
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

    while True:
        clusters = [[] for c in centroids]
        for p in points:
            index = 0
	    # calculate the distance from that point to all centroids
            smallest_distance = float('inf')
            for i in range(len(centroids)):
                distance = getDistance(p, centroids[i].centroid)
                if distance < smallest_distance:
                    smallest_distance = distance
                    index = i
	    # the point belongs to the centroid of least distance to it
            clusters[index].append(p)
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

