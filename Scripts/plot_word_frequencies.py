"""
This function will plot the frequencies of input
words in Medline titles over a specified time period. It only displays words that were used 
in a title at least 100 times during the study period (1850-Aug 2014). The first two inputs
are years and all following inputs are words. Sample call (rat,mouse,chicken from 1980 to 2010,inclusive):
plot_word_frequencies.py 1980 2010 rat mouse chicken
"""

import pandas as pd
import pickle
import matplotlib.pyplot as plt
import sys


if __name__ == "__main__":
   	y1 = sys.argv[1] #starting year
   	y2 = sys.argv[2] #ending year
   	inputs = sys.argv[3:] #words to query
   	     
	#load word frequencies
	WP = pickle.load(open('/Users/jschwarz/Code/Python/Projects/Processed data/Pubmed/WordFreq','r'))
	WP_index = WP.index
	WP_years = WP.columns

	#Check that all input words are in the dataset, only include those that are
	ind = [] #indices of input words present in dataset
	for i in inputs:
		if i in WP_index:
			ind.append(i)
	#Check that both years are in dataset
	if y1 not in WP_years: y1 = WP_years[0]
	if y2 not in WP_years: y2 = WP_years[-1]
	
	WP.ix[ind,y1:y2].T.plot() #plot the subset of the data
	plt.show()