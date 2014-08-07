"""This script will load word_counts and bigram_counts from pickle files across years and
place them in data frames, then pickle them
"""

import pandas as pd
from os import listdir
import numpy as np
import pickle

processed_data_dir = '../../Processed data/Pubmed/' #location of processed data

#Load word frequencies
file_list = listdir(processed_data_dir + 'Words/')

#make sure files are in the correct order. File name format is words_YYYY
file_list.sort()

word_dict = {}

for f in file_list:
	print('Loading file {fn}').format(fn=f)
	p_file = open(processed_data_dir + 'Words/' + f,'r')
	word_dict[f[-4:]] = pickle.load(p_file)
	p_file.close()

print('Building WordFreq DataFrame')	
WordFreq = pd.DataFrame(word_dict)
row_sums = WordFreq.apply(lambda x: np.sum(x),axis=1)
WordFreq_slim = WordFreq[row_sums>100]
WordFreq_slim.fillna(0,inplace=True) #replace NaN with 0

WordProp = WordFreq_slim.apply(lambda x: x/np.sum(x))

p_file = open(processed_data_dir + 'WordFreq','w')
pickle.dump(WordFreq_slim,p_file)
p_file.close()

p_file = open(processed_data_dir + 'WordProp','w')
pickle.dump(WordProp,p_file)
p_file.close()
"""
#Load bigram frequences
file_list = listdir(processed_data_dir + 'Bigrams/')

#make sure files are in the correct order. File name format is words_YYYY
file_list.sort()

bigram_dict = {}

for f in file_list:
	print('Loading file {fn}').format(fn=f)
	p_file = open(processed_data_dir + 'Bigrams/' + f,'r')
	word_dict[f[-4:]] = pickle.load(p_file)
	p_file.close()

print('Building BigramFreq DataFrame')	
BigramFreq = pd.DataFrame(bigram_dict)
row_sums = BigramFreq.apply(lambda x: np.sum(x),axis=1)
BigramFreq_slim = BigramFreq[row_sums>100]
BigramFreq_slim.fillna(0,inplace=True) #replace NaN with 0


p_file = open(processed_data_dir + 'BigramFreq','w')
pickle.dump(BigramFreq_slim,p_file)
p_file.close()
"""