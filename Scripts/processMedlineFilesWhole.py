"""
This script will process .txt files in the specified directory and save within subdirectories
the word and bigram frequencies for each file
"""

from os import listdir
from pyspark import SparkContext, SparkConf
import pandas as pd

sc = SparkContext()

def getMedlineTitleWord(data_directory):
	"""
	pulls titles from directory of *.txt files containing Medline entries and finds the
	frequency of words and bigrams within those titles for each file
	"""
	def removePunctuation(x):
		"""This function removes various punctuation and lowercases all strings"""
		
		return x.replace('.','').replace(',','').replace('?','').replace('-',' ').replace('   ',' ') \
		.replace('   ',' ').replace('"','').replace("'",'').replace('(','').replace(')','').replace(':','') \
        .replace('[','').replace(']','').replace('\\','').replace('/','').replace('!','').replace('/n','').lower()
	
	"""
	Extract titles from Medline entries:
	1.Make RDD from text files;	2. Split strings on PMID flag; 3. Pull text following TI flag;
	4. Pull text up to next field flag; 5. Remove punctuation, extra spaces, make lower case
	titles is an RDD consisting of a list of 2-tuples, the first element is the file name and
	the second element is each title
	"""
	titles = sc.textFile(text_file) \
	.flatMapValues(lambda x: x.split("PMID- ")) \
	.mapValues(lambda x: x[x.find('TI ')+5:]) \
	.mapValues(lambda x: x[:(x.find('  - ')-3)]) \
	.mapValues(removePunctuation)
	
	"""
	Get word frequencies from titles:
	1. Split titles into single words; 2. make each word a 2-element tuple, with the word in 
	first position and integer 1 in second; 3. reduce by finding all unique tuple[0] (word) values;
	4. Collect RDD as dictionary
	"""
	title_string = titles.reduceByKey(lambda x,y:x+1) \
	.combineAsMap()
	
	title_word_freq = titles.flatMap(lambda x: x.split()) \
	.map(lambda x: (x,1)) \
	.reduceByKey(lambda x,y:x+y) \
	.collectAsMap()
	
	"""
	Get title word bigram frequencies: 
	1. Break title string into individual titles; 2. For each title, return all word pairings 
	as a 3-tuple: word1,word2,1; 3. reduce by finding all unique tuple[0:1] (bigram) values; 
	4. Collect RDD as dictionary
	"""
	title_bigram_freq = titles.reduceByKey(lambda x,y: x + ',' + y) \
	.mapValues(lambda x: x.split(',')) \
	
	.flatMap(lambda x: [((x[i],x[i+1]),1) for i in range(0,len(x)-1)]) \
	.reduceByKey(lambda x,y:x+y) \
	.collectAsMap()
	
	return title_word_freq, title_bigram_freq


data_directory = '../../Data/Pubmed/medline/' #Folder containing medline data files in .txt format

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