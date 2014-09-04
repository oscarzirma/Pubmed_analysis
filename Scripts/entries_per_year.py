"""
This script will use the Entrez programming utility to find the number of Pubmed entries by year
"""
from Bio import Entrez
from Bio import Medline

start_year = 1850 #first year from which to pull medline entries (inclusive)
end_year =  2014 #last year from which to pull medline entries (inclusive)

data_directory = '../../Data/Pubmed/medline/'

Entrez.email = "schwarzj@stanford.edu" #email for entrez tracking

#non_year_search_term = "(Smith[AUTH]) AND (" #search term preceding the year (example)
non_year_search_term = "" #for this data set, we only use the year
num_entries = {}
for y in range(start_year,end_year+1):
	print("Year: %i" % y)
	search_results = Entrez.read(Entrez.esearch(db="pubmed",term=non_year_search_term + str(y) + "[PDAT])",usehistory="y")) #pubmed query
	count = int(search_results["Count"]) #number of query results
	num_entries[y] = count

f = open(data_directory + 'numberEntries','w')
f.writelines(str(num_entries))
f.close()