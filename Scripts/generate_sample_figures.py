"""
This script generates, annotates, and saves a set of sample figures from the medline title
word frequency analysis
"""

import pandas as pd
import pickle
import matplotlib.pyplot as plt
import sys

def checkInputs(inputs,y1,y2):
	#Check that all input words are in the dataset, only include those that are
	ind = [] #indices of input words present in dataset
	for i in inputs:
		if i in WP_index:
			ind.append(i)
	#Check that both years are in dataset
	if y1 not in WP_years: y1 = WP_years[0]
	if y2 not in WP_years: y2 = WP_years[-1]
	
	return ind,y1,y2

	 
#load word frequencies
WP = pickle.load(open('../../Processed data/Pubmed/WordFreq','r'))
WP_index = WP.index
WP_years = WP.columns


inputs = ['influenza']#words to query
year1 = '1870'
year2 = '2013'

ind,y1,y2 = checkInputs(inputs,year1,year2

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(WP.ix[ind,y1:y2].T*100) #plot the subset of the data
ax.xlabel('Year')
ax.ylabel('Percent of title words')
ax.savefig('../Figures/' + inputs[0] + '.eps')