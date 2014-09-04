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
	
	def journal_title_key_val(x):
		"""This function finds title and journal id fields in medline entry"""
		ti1 = x.find('TI ')+5
		ti2 = x[ti1:].find('  - ')-3
		ji1 = x.find('JID ')+6
		ji2 = x[ji1:].find('- ')-5
		return (x[ji1:ji1+ji2],removePunctuation(x[ti1:ti1+ti2]))
	
	"""Separate file into seperate medline entries"""
	entries = sc.textFile(text_file) \
	.glom() \
	.map(lambda x: " ".join(x)) \
	.flatMap(lambda x: x.split("PMID- "))
	
	
	"""Join titles from same journal"""
	titles_by_journal = entries.map(journal_title_key_val) \
	.reduceByKey(lambda x,y:x+y) \
	.mapValues(lambda x:x.split()).collectAsMap()
	
	titles_by_journal.pop('')
	"""For each journal, get word counts"""
	journal_word_freq={}
	for i in titles_by_journal:
		print i
		journal_word_freq[i] = sc.parallelize(journal_word_freq[i]) \
		.flatMap(lambda x: x.split()) \
		.map(lambda x: (x,1)) \
		.reduceByKey(lambda x,y:x+y) \
		.collectAsMap()
		
	return journal_word_freq

"""
data_directory = '../../Data/Pubmed/medline/' #Folder containing medline data files in .txt format
processed_data = '../../Processed data/Pubmed/' #location to save processed data
"""
data_directory = '/Users/jschwarz/Code/Python/Projects/Data/Pubmed/test/' #Folder containing medline data files in .txt format
processed_data = '/Users/jschwarz/Code/Python/Projects/Processed data/test/' #location to save processed data
file_list = listdir(data_directory) #get list of all files in data_directory

#initialize frequency dictionaries
word_dict = {}

for f in file_list:
	if f[-3:] == 'txt': #if file is txt file	
		print('Current file: {fn}. word_list length: {lwf}'.format(fn=f,lwf=len(word_dict))) #Print current file name and dict size
		word_freq = getMedlineTitleWord(data_directory + f) #Run spark script to process files and get word,bigram frequencies
		word_dict[f[-9:-4]] = word_freq #add word frequencies from current file to word_dict keyed to subset of current file name
		bigram_dict[f[-9:-4]] = bigram_freq #add bigram frequencies from current file to word_dict keyed to subset of current file name
		print(word_dict.keys()) #list all of the current keys
		output_file = open(processed_data+'/words'+f[-9:-4],'w')
		pickle.dump(word_freq,output_file)
		output_file.close()

