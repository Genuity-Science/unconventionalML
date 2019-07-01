#####################################################################################################
##### PMC Crawling for WuXi NextCODE AI lab
##### Code: Sweta Bajaj
#####################################################################################################

# Load the Libraries required
options(stringsAsFactors = FALSE)
library("rentrez")
library("hmisc")

# Check the Database Summary
entrez_db_searchable("pmc")
entrez_db_summary("pmc")

#####################################################################################################
# Load the Data for the gene list inputs
setwd("~/OneDrive - NextCODE Health/Work_Experiments/quantum/NLP_tool_PMC/")
input_gene_list <- read.csv("LumA_LumB_input_forR.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#####################################################################################################
# Set Parameters
# Add the databse name from which the query terms are to be searched
db_name <- "pmc"
# Add the second term with the gene : eg. cancer or breast canver
second_search_term <- "breast cancer"
#####################################################################################################
# Query the database
paper_names_list <- list()
gene_hits <- list()
for (i in 1:nrow(input_gene_list)){
  # Create the Query Term
  query <- paste0(input_gene_list$V1[i],"[gene] and ", second_search_term )
  
  # Search for the query term and Cancer term
  r_search <- entrez_search(db=db_name, term=query , retmax=50)

   # if the query term has any results >0
  if(length(r_search$ids) > 0) {

    # Get the paper Names
    article_summary <- entrez_summary(db=db_name, id=r_search$ids)
    article_journalnames <- extract_from_esummary(article_summary, "fulljournalname")
    article_title <- extract_from_esummary(article_summary, "title")
    
    article_names_table <- do.call(rbind, Map(data.frame,"Journal_Name" = article_journalnames, "Title" = article_title))
    
  paper_names_list[[i]] <- article_names_table
  } else {
     paper_names_list[[i]] <- "No hits for this Gene"
  }
  
  names(paper_names_list)[i] <- input_gene_list$V1[i]
  gene_hits[[i]] = entrez_search(db=db_name, term=query)$count
}

# Create a dataframe to Enter the Gene and Hit counts 
summary_df <- do.call(rbind, Map(data.frame, 'Gene_Name'= input_gene_list$V1,'Hits' = gene_hits))
write.csv(summary_df, paste0(db_name, "_", second_search_term, "_Gene_Link_Counts.csv"), row.names = FALSE, quote = FALSE)
save(paper_names_list, file = paste0(db_name, "_", second_search_term, "_Gene_Links.RData"))
#####################################################################################################

#####################################################################################################
# Create the bins for Circos Plots
#Load the cancer link files
Cancer_Gene_Link_Counts <- read.csv("pmc_cancer_Gene_Link_Counts_2.csv", sep = ",")
breastCanvcr_Gene_Link_Counts <- read.csv("pmc_breast cancer_Gene_Link_Counts_2.csv", sep = ",")
cancer_hits_bins <- cbind(Cancer_Gene_Link_Counts, Hmisc::cut2( Cancer_Gene_Link_Counts$Hits, g=5))
View(cancer_hits_bins)
breast_cancer_hits_bins <- cbind(breastCanvcr_Gene_Link_Counts, Hmisc::cut2( breastCanvcr_Gene_Link_Counts$Hits, g=5))
View(breast_cancer_hits_bins)
#####################################################################################################
