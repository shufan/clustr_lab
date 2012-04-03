import sys
import csv
import getopt
import random

bases = ['A', 'C', 'G', 'T']

def usage():
   print '$> python dnagen.py <required args>\n' + \
      '\t-c <#>\t\tNumber of clusters to generate\n' + \
      '\t-p <#>\t\tNumber of strands per cluster\n' + \
      '\t-o <file>\tFilename for the output of the raw data\n' + \
      '\t-v [#]\t\tLength of the DNA strand\n'

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
def maxSimilarity(centroid, centroidArr, dnaLength):
   maxSim = 0
   for pair in centroidArr:
      sim = similarity(centroid, pair, dnaLength)
      if sim > maxSim:
	 maxSim = sim
   return maxSim

#creates a modified DNA strand based on a given centroid
def createStrand(centroid, dnaLength):
   newStrand = list(centroid)
   numDiffs = int(dnaLength * .25)
   for i in range(numDiffs):
      changeIndex = int(random.random()*dnaLength)
      newStrand[changeIndex] = random.choice(bases)
   return newStrand

def check(centroidArr, strandsArr):
   for strand in strandsArr:
      #print strand
      for centroid in centroidArr:
	 count = 0
	 for i in range(len(centroid)):
	    if centroid[i] == strand[i]:
	       count += 1
	 print count, 
      print "\n",

def handleArgs(args):
   #set up return values
   numClusters = -1
   numPerCluster = -1
   output = None
   dnaLength = 10

   try:
      optlist, args = getopt.getopt(args[1:], 'c:p:v:o:')
   except getopt.GetoptError, err:
      print str(err)
      usage()
      sys.exit(2)
   for key, val in optlist:
      #first, the required arguments
      if   key == '-c':
         numClusters = int(val)
      elif key == '-p':
         numPerCluster = int(val)
      elif key == '-o':
         output = val
      #now, the optional argument
      elif key == '-v':
	 dnaLength = int(val)

   #check required arguments were inputted
   if numClusters < 0 or numPerCluster < 0 or dnaLength < 1 or output is None:
      usage()
      sys.exit()

   return (numClusters, numPerCluster, output, dnaLength)

def main():
   #start by reading the command line
   numClusters, \
   numPerCluster, \
   output, \
   dnaLength = handleArgs(sys.argv)

   writer = csv.writer(open(output, "w"))

   #step 1: create centroid DNA strands
   centroidArr = []
   maxAllowedSim = int(dnaLength * .75)
   for i in range(numClusters):
      centroid = createCentroid(dnaLength)
      tryCount = 0
      #is it different enough from other centroids in centroidArr?
      similarity = maxSimilarity(centroid, centroidArr, dnaLength)
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
	 similarity = maxSimilarity(centroid, centroidArr, dnaLength)
	 tryCount += 1
      centroidArr.append(bestCentroid)

   #step 2: create groups of numPerCluster strands each based on the centroids
   strandsArr = []
   for centroid in centroidArr:
      for i in range(numPerCluster):
	 newStrand = "".join(createStrand(centroid, dnaLength))
	 while (newStrand in strandsArr):
	     newStrand = "".join(createStrand(centroid, dnaLength))
	 strandsArr.append(newStrand)

   for strand in strandsArr:
      writer.writerow([strand])

#   check(centroidArr, strandsArr)

if __name__ == "__main__":
    main()
