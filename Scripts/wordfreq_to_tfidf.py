"""
This script loads WordFreq and WordProp dataframes and calculates tf*idf for each word
"""

import pandas as pd
import pickle
import numpy as np

processed_data_dir = '../../Processed data/Pubmed/' #location of processed data

print('Load pickles')
WordFreq = pickle.load(open(processed_data_dir + 'WordFreq','r'))
WordProp = pickle.load(open(processed_data_dir + 'WordProp','r'))

print('Calculate WordTfIdf')
WordPresence = WordFreq>0 #dataframe of word presence and absence in each year

document_count = WordPresence.apply(np.sum,axis=1) #number of years in which each word is present

num_words,num_years = WordPresence.shape

inverse_document_freq = np.log10((num_years+.1-.1)/document_count) #calculate the inverse document frequency for each word

WordTfIdf = WordProp.apply(lambda x: np.asarray(x) * np.asarray(inverse_document_freq))

print('Save WordTfIdf')
p_file = open(processed_data_dir + 'WordTfIdf','w')
pickle.dump(WordTfIdf,p_file)
p_file.close()