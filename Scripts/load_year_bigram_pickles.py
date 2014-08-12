"""This script will load bigram_counts from pickle files across years and
place them in data frames, then pickle them
"""

import pandas as pd
from os import listdir
import numpy as np
import pickle

processed_data_dir = '../../Processed data/Pubmed/' #location of processed data


#Load bigram frequences
file_list = listdir(processed_data_dir + 'Bigrams/')

#make sure files are in the correct order. File name format is words_YYYY
file_list.sort()

BigramFreq = pd.DataFrame([])

for f in file_list:
	print('Loading file {fn}').format(fn=f)
	p_file = open(processed_data_dir + 'Bigrams/' + f,'r') #load each bigram pickle
	word_dict={f[-4:] : pickle.load(p_file)}
	p_file.close()
	BigramFreq = pd.concat([BigramFreq,pd.DataFrame(word_dict)] #add each bigram pickle to data frame
	del word_dict

#print('Building BigramFreq DataFrame')	
#BigramFreq = pd.DataFrame(bigram_dict)
#row_sums = BigramFreq.apply(lambda x: np.sum(x),axis=1)
#BigramFreq_slim = BigramFreq[row_sums>100]
#BigramFreq_slim.fillna(0,inplace=True) #replace NaN with 0


p_file = open(processed_data_dir + 'BigramFreq','w')
pickle.dump(BigramFreq_slim,p_file)
p_file.close()