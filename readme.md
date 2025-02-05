Analysis of Medline titles over time
==========
What can we learn about science by studying title word frequencies over time?
------------

The goal of this project is to think about the information provided by the titles of scientific papers. 
My hope is that we can learn about trends in science and scientific thinking.
To see reports generated by these scripts see my [webpage](http://web.stanford.edu/~schwarzj/)

##Python scripts, in order of pipeline
**pull_pubmed_citations** - pull Medline entries by year and query  
**process_medline_files** - process Medline entries into title word and bigram frequencies by year  
**load_year_freq_pickles** - loads word frequency pickle dicts into pandas dataframe  
**plot_word_frequencies** - loads word frequency dataframe and plots specified words in specified years  
**plot_word_frequencies_and_save** - loads word frequency dataframe and saves .eps plot of specified words in specified years  
**word_freq_chi** - does a chi-squared analysis on word frequencies as a first pass look at trending words (incomplete)  