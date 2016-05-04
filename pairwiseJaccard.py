#!usr/bin/env python

#pairwiseJaccard.py
#Author: Mark Hartman
#Date created: 3/31/2016
#Last modified: 5/4/2016

#Given a spreadsheet of counts for various genus (or other OTU level) for each sample
# Format as shown:
#Genus,MH01-A,MH95HP,MH10-B, ...
#Streptococcus,38,79115,6225, ...
#Propionibacterium,26,3,36045, ...
#Staphylococcus,5,3,6028, ...
#...
# (Note: This type of file is available from the BaseSpace 16S Metagenomics App.)
#
#
#Return a pairwise Jaccard calculation
# Format as shown:
#[blank cell], MH01-A, MH95HP, MH10-B, ...
#MH01-A, 1, 0.111111111, 0.25, ...
#MH95HP, 0.111111111, 1, 0.25, ...
#MH10-B, 0.25, 0.25, 1, ...
#...

INPUT_FILENAME = 'Genus_Level_Aggregate_Counts.csv'
OUTPUT_FILENAME = 'PairwiseJaccardResults.csv'

def main():
	readDataOutput = readData(INPUT_FILENAME)
	listOfSamples = readDataOutput[0]
	listOfGenus = readDataOutput[1]
	aggregateCounts = readDataOutput[2]
	dictOfSampleNamesTopGenus = makeDictOfTopGenus(aggregateCounts, listOfSamples, listOfGenus)
	listOfJaccards = getJaccards(dictOfSampleNamesTopGenus) 
	writeJaccards(listOfSamples, listOfJaccards, OUTPUT_FILENAME) 

def readData(INPUT_FILENAME):
	#read in sample names, genus names, and aggregate counts
	with open(INPUT_FILENAME, 'r') as filer:
		listOfSamples = filer.readline().strip().split(',')[1:]

		listOfGenus = []
		aggregateCounts = []
		for line in filer:
			currentLine = line.strip().split(',')
			listOfGenus.append(currentLine[0])
			aggregateCounts.append(currentLine[1:])
	aggregateCounts = map(list, zip(*aggregateCounts)) #transpose the aggregate counts
	aggregateCounts = [map(int, x) for x in aggregateCounts] #convert aggregate counts to int
	return (listOfSamples, listOfGenus, aggregateCounts)

#convert aggregate counts to presence/absence calls for each genus
#make a list of sets, each set contains the ten most abundent genus in each sample
def makeDictOfTopGenus(aggregateCounts, listOfSamples, listOfGenus):
	listOfMostAbundantGenus = []
	for samplewiseAggregateCounts in aggregateCounts:
		listOfMaxValueIndexes = []
		mostAbundantGenus = set() #initializes the set that will contain the 10 most abundent genus
		for x in range(10):
	 		currIndex = samplewiseAggregateCounts.index(max(samplewiseAggregateCounts))
	 		samplewiseAggregateCounts[currIndex] = 0
	 		mostAbundantGenus.add(listOfGenus[currIndex])
	 	listOfMostAbundantGenus.append(mostAbundantGenus)
	#zip these sets into a dictionary {key = sample name: value = set of most abundant genus}
	dictOfSampleNamesTopGenus = dict(zip(listOfSamples,listOfMostAbundantGenus))
	return dictOfSampleNamesTopGenus

#calculate the Jaccard similarity for each pair of samples
def getJaccards(dictOfSampleNamesTopGenus):
	listOfJaccards = []
	for firstEntry in dictOfSampleNamesTopGenus:
		for secondEntry in dictOfSampleNamesTopGenus:
			setA = dictOfSampleNamesTopGenus[firstEntry]
			setB = dictOfSampleNamesTopGenus[secondEntry]
			jaccardSimilarity = (float(len(setA&setB))/float(len(setA|setB)))
			listOfJaccards.append(jaccardSimilarity)
	listOfJaccards = map(str, listOfJaccards)
	return listOfJaccards

#writes the Jaccard similarities for each pair of samples to a csv
def writeJaccards(listOfSamples, listOfJaccards, OUTPUT_FILENAME):
	with open(OUTPUT_FILENAME, 'w') as filew:
		filew.write(',')
		filew.write(','.join(listOfSamples))
		filew.write('\n')
		numOfSamples = len(listOfSamples)
		for i in range(numOfSamples):
			filew.write(listOfSamples[i] + ',')
			filew.write(','.join(listOfJaccards[i*numOfSamples:(i+1)*numOfSamples]))
			filew.write('\n')
	print 'Jaccard\'s printed to', OUTPUT_FILENAME

if __name__ == '__main__':
        main()

