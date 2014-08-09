"""
This script will process .txt files in the specified directory and save within subdirectories
the word and bigram frequencies for each file
"""

from os import listdir
from pyspark import SparkContext, SparkConf
import pandas as pd
import pickle

sc = SparkContext()

def getMedlineTitleWord(text_file):
	"""
	pulls titles from .txt file containing Medline entries and finds the
	frequency of words and bigrams within those titles
	"""
	def removePunctuation(x):
		"""This function removes various punctuation and lowercases all strings"""
		
		return x.replace('.','').replace(',','').replace('?','').replace('-',' ').replace('   ',' ') \
		.replace('   ',' ').replace('"','').replace("'",'').replace('(','').replace(')','').replace(':','') \
        .replace('[','').replace(']','').replace('\\','').replace('/','').replace('!','').replace('~','') \
        .replace('{','').replace('}','').replace('#','').lower()
	
	"""
	Extract titles from Medline entries:
	1.Make RDD from text files; 2. Glom files together; 3. Join glommed entries into string; 
	4. Split string on PMID flag; 5. Pull text following TI flag; 6. Pull text up to next field flag; 
	7. Remove punctuation, extra spaces, make lower case
	"""
	titles = sc.textFile(text_file) \
	.glom() \
	.map(lambda x: " ".join(x)) \
	.flatMap(lambda x: x.split("PMID- ")) \
	.map(lambda x: x[x.find('TI ')+5:]) \
	.map(lambda x: x[:(x.find('  - ')-3)]) \
	.map(removePunctuation)
	
	"""
	Get word frequencies from titles:
	1. Split titles into single words; 2. make each word a 2-element tuple, with the word in 
	first position and integer 1 in second; 3. reduce by finding all unique tuple[0] (word) values;
	4. Collect RDD as dictionary
	"""
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
	title_bigram_freq = titles.map(lambda x:x.split()) \
	.flatMap(lambda x: [((x[i],x[i+1]),1) for i in range(0,len(x)-1)]) \
	.reduceByKey(lambda x,y:x+y) \
	.collectAsMap()
	
	return title_word_freq, title_bigram_freq


data_directory = '../../Data/Pubmed/medline/' #Folder containing medline data files in .txt format
processed_data = '../../Processed data/Pubmed/' #location to save processed data

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
		output_file = open(processed_data+'/words'+f[-9:-4],'w')
		pickle.dump(word_freq,output_file)
		output_file.close()
		output_file = open(processed_data+'/bigrams'+f[-9:-4],'w')
		pickle.dump(bigram_freq,output_file)
		output_file.close()		
