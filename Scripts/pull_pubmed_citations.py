"""
This script will use the Entrez programming utility to pull PubMed searches by year.
For each year, it creates a .txt file containing the Medline entries. The script will evaluate
go through the ./medline directory and check for complete years based on medline_log.txt,
only downloading those years which are incomplete or nonexistant. See Entrez e-utility documentation
and Biopython Entrez methods for further information on query structure.
"""
from Bio import Entrez
from Bio import Medline

start_year = 1850 #first year from which to pull medline entries (inclusive)
end_year =  1940 #last year from which to pull medline entries (inclusive)

data_directory = '../../Data/Pubmed/medline/'

batch_size = 2000;#record to pull in each iteration
Entrez.email = "schwarzj@stanford.edu" #email for entre tracking
#non_year_search_term = "(Smith[AUTH]) AND (" #search term preceding the year (example)
non_year_search_term = "" #for this data set, we only use the year

#Open log file and generate a list to track completed years
log_file = open(data_directory + 'medline_log','r')
completed_years = log_file.readlines()
log_file.close()
log_file = open(data_directory + 'medline_log','a')
incomplete = True

for y in range(start_year,end_year+1):
	print("Year: %i" % y)
	#Find if year has already been pulled
	for cy in completed_years:
		if str(y) == cy.rstrip():
			incomplete=False
			print("Year already complete")
			break
		else:
			incomplete=True
	
	#If year has not been pulled, do it
	if incomplete:
		file_name = data_directory + 'medline_data_' + str(y) + '.txt' #destination file for data
		search_results = Entrez.read(Entrez.esearch(db="pubmed",term=non_year_search_term + str(y) + "[PDAT])",usehistory="y")) #pubmed query
		count = int(search_results["Count"]) #number of query results
		print("Found %i results" % count)
		out_handle = open(file_name,'w') #open data destination file
		for start in range(0,count,batch_size):#for each batch, fetch medline query results
			end = min(count,start+batch_size)
			print("Going to download record %i to %i out of %i" % (start+1, end, count))
			fetch_handle = Entrez.efetch(db="pubmed",
                                 rettype="medline", retmode="csv",
                                 retstart=start, retmax=batch_size,
                                 webenv=search_results["WebEnv"],
                                 query_key=search_results["QueryKey"]) #Fetch command
			data = fetch_handle.read() #read handle
			fetch_handle.close() #close handle
			out_handle.write(data) #write data to data destination file
		out_handle.close() #close data destination file
		log_file.writelines(str(y) + "\n") #log this year as completed

log_file.close() #close log file