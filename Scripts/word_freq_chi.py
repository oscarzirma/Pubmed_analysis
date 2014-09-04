"""
This script will process WordFreq and WordProp to produce expected counts for each word 
based on overall proportion and words per year. It will then calculate a Chi-square statistic
by squaring the difference between expected and observed and dividing it by expected value
for each word in each year
"""

import pandas as pd
import pickle
import numpy as np

WF = pickle.load(open('/Users/jschwarz/Code/Python/Projects/Processed data/Pubmed/WordFreq','r'))
WP = pickle.load(open('/Users/jschwarz/Code/Python/Projects/Processed data/Pubmed/WordProp','r'))

words = WF.index #row labels (indices), which are the words
years = WF.columns #column labels, which are the years

total_words = WF.apply(sum,axis=1) #for each word, how many total instances
total_words_prop = total_words/sum(total_words) #for each word, what proportion is it of all words?
words_per_year = WF.apply(sum) #how many total words in each year

num_words = len(total_words_prop) #number of different words in the dataset
num_years = len(word_per_year) #number of years in the dataset

WFexp = np.ndarray((num_words,num_years)) #initialize matrix

#Build the expected matrix by multiplying, for each word in each year, its overall expected proportion times the total words in that year (vectorize this code?)
for i in np.arange(0,num_words):
    for j in np.arange(0,num_years):
        WFexp[i,j] = word_per_year[j]*total_words_prop[i]
        
WFexpdf = pd.DataFrame(WFexp,index=words,columns=years) #Create a dataframe from matrix of expected word count

WFexp_diff = WF-WFexpdf #Find the difference between observed and expected counts

WFexp_diff_sq = WFexp_diff.apply(np.square) #square the difference between observed and expected counts

WF_chi = WFexp_diff_sq/WFexpdf #divide the squared difference by the expected counts

word_chi_val = WF_chi.apply(np.sum,axis=1) #for each word, find the total chi statistic


##Save files