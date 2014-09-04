"""
This script will process .txt files in the specified directory and save within subdirectories
an array of tuples: (journalID,title)
"""

from os import listdir
from pyspark import SparkContext, SparkConf
import pandas as pd
import pickle

sc = SparkContext()

def getMedlineTitleJournalID(text_file):
	"""
	pulls titles and journalIDs from .txt file containing Medline entries
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
		jti1 = x.find('JT ')+5
		jti2 = x[jti1:].find(' - ')-4
		return (x[ji1:ji1+ji2],x[jti1:jti1+jti2],removePunctuation(x[ti1:ti1+ti2]))
	
	"""Separate file into seperate medline entries"""
	entries = sc.textFile(text_file) \
	.glom() \
	.map(lambda x: " ".join(x)) \
	.flatMap(lambda x: x.split("PMID- "))
	
	
	"""Join titles from same journal"""
	titles_by_journal = entries.map(journal_title_key_val)
	
	titles_by_journal.pop('')
		
	return journal_word_freq



"""
data_directory = '../../Data/Pubmed/medline/' #Folder containing medline data files in .txt format
processed_data = '../../Processed data/Pubmed/' #location to save processed data
"""
data_directory = '/Users/jschwarz/Code/Python/Projects/Data/Pubmed/test/' #Folder containing medline data files in .txt format
processed_data = '/Users/jschwarz/Code/Python/Projects/Processed data/Journals/' #location to save processed data
file_list_full = listdir(data_directory) #get list of all files in data_directory
file_list = file_list_full[148:]

for f in file_list:
	if f[-3:] == 'txt': #if file is txt file	
		print('Current file: {fn}.'.format(fn=f)) #Print current file name and dict size
		journalID_title = (getMedlineTitleJournalID(data_directory + f)) #Run spark script to process files and get word,bigram frequencies
		print('Writing file')
		output_file = open(processed_data+'journals'+f[-9:-4],'w')
		pickle.dump(journalID_title,output_file)
		output_file.close()

