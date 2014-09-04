jid_t <- read.csv(file = "~/Code/Python//Projects/Processed data/Pubmed/journalID_title_1990.csv",stringsAsFactors=FALSE)
impact <- read.csv(file = "~/Code/Python//Projects/Data/Pubmed//SJR journal rankings and impact.tsv",sep='\t',stringsAsFactors=FALSE)

tempcol <- matrix(c(rep.int(NA,nrow(jid_t))),nrow = nrow(jid_t),ncol=1)
tempdf <- data.frame(tempcol)
colnames(tempdf)<-'SJR'
journals <- cbind(jid_t,tempdf)

colnames(tempdf)<-'H.Index'
journals <- cbind(journals,tempdf)