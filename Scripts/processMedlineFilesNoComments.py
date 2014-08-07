"""
This script will process .txt files in the specified directory and save within subdirectories
the word and bigram frequencies for each file
"""

from os import listdir
from pyspark import SparkContext, SparkConf
import pandas as pd

sc = SparkContext()

def getMedlineTitleWord(text_file):

	def removePunctuation(x):
		
		return x.replace('.','').replace(',','').replace('?','').replace('-',' ').replace('   ',' ') \
		.replace('   ',' ').replace('"','').replace("'",'').replace('(','').replace(')','').replace(':','') \
        .replace('[','').replace(']','').replace('\\','').replace('/','').replace('!','').lower()
	
	titles = sc.textFile(text_file) \
	.glom() \
	.map(lambda x: " ".join(x)) \
	.flatMap(lambda x: x.split("PMID- ")) \
	.map(lambda x: x[x.find('TI ')+5:]) \
	.map(lambda x: x[:(x.find('  - ')-3)]) \
	.map(removePunctuation)
	
	title_word_freq = titles.flatMap(lambda x: x.split()) \
	.map(lambda x: (x,1)) \
	.reduceByKey(lambda x,y:x+y) \
	.collectAsMap()
	
	title_bigram_freq = titles.map(lambda x:x.split()) \
	.flatMap(lambda x: [((x[i],x[i+1]),1) for i in range(0,len(x)-1)]) \
	.reduceByKey(lambda x,y:x+y) \
	.collectAsMap()
	
	return title_word_freq, title_bigram_freq


data_directory = '../../Data/Pubmed/test/' #Folder containing medline data files in .txt format

file_list = listdir(data_directory) #get list of all files in data_directory

#initialize frequency dictionaries
word_dict = {}
bigram_dict = {}

for f in file_list:
	if f[-3:] == 'txt': #if file is txt file
		print('Current file: {fn}. word_list length: {lwf}'.format(fn=f,lwf=len(word_dict))) #Print current file name and dict size
		word_freq, bigram_freq = getMedlineTitleWord(data_directory + f) #Run spark script to process files and get word,bigram frequencies
		word_dict[f[-9:-4]] = word_freq #add word frequencies from current file to word_dict keyed to subset of current file name
		bigram_dict[f[-9:-4]] = bigram_freq #add bigram frequencies from current file to word_dict keyed to subset of current file name
		print(word_dict.keys()) #list all of the current keys
WordFreq = pd.DataFrame.from_dict(word_dict) #generate data frame from word frequency dict
BigramFreq = pd.DataFrame.from_dict(bigram_dict) #generate data frame from bigram frequency dict

WordFreq.fillna(0,inplace=True) #replace NaN with 0
WordFreq.to_pickle(data_directory + 'WordFreq_pickle') #save dataframe in pickle format

BigramFreq.fillna(0,inplace=True) #replace NaN with 0
BigramFreq.to_pickle(data_directory + 'BigramFreq_pickle') #save dataframe in pickle format