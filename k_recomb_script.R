library(RMySQL)
library(arules)
library(arulesViz)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggplot2)
library(data.tree)
library(igraph)
library(viridis)
library(qgraph)
library(FSA)
library(rcompanion)

#********************************************************************************************************************************************************************************************************
# Code to investigate knowledge recombination novelty in the invention stages of the innovation process and its relation with collaboration breadth. It builds on IPC codes combinations in patent. 
# STEP 1:
#     1A: Get patents data from DB
#     1B: Inspect data for missing values, zeros, extreme outliers or erroneous characters
# STEP 2: Run association rule mining     
# STEP 3: Construct novelty measure
# STEP 4: Clean the applicants column
# STEP 5: Create the collaboration variable
# STEP 6: Run ANOVA
# STEP 7: Visualization
#********************************************************************************************************************************************************************************************************
options(scipen=999) #disable scientific notation

#===STEP 1 (get data and inspect) =============================================================================================================================================================
#1A. Get data from DB
psswd <- .rs.askForPassword("Database password:") #activate a pop up asking for the password, which is then stored in this variable
mydb <- RMySQL::dbConnect(MySQL(), user = 'root', password = psswd, dbname = 'CPS_project_patents') # connect with MySQL database
query <- function(...) dbGetQuery(mydb, ...) #---Function to make queries faster and easier
patents <- query("SELECT id, pub_date, country, applicants, inventors FROM patents WHERE patents.CPS_group = 'Robotics and Autonomous Systems';")
patents_ipc_code <- query("SELECT patent_id AS 'id', ipc_class, ipc_subclass, ipc_maingroup FROM patents_ipc_codes")
data <- merge(patents, patents_ipc_code, by = "id") # merge data

#1B. Inspect and clean data
# Most of the inspection for erroneous characters and missing IPC codes was done during data insertion into db
data <- data[!is.na(data$applicants) & !is.na(data$inventors), ] #filter out patents with missing applicants and/or inventors

data <- data %>% 
  filter(str_detect(applicants, "__", negate = T)) #remove patents with patterns like "__" in the applicants (this is likely to be text in Mandarin, Japanese, Korean, etc.)

data <- data %>% 
  filter(str_detect(inventors, "__", negate = T)) #remove patents with patterns like "__" in the applicants (this is likely to be text in Mnadarin, Japanese, Korean, etc.)
#===End STEP 1=================================================================================================================================================================================




#===STEP 2 (Run association rule mining) =============================================================================================================================================================
# Format and pre-process the data:
data2019 <- data[substr(data$pub_date,7,10) == '2019',c(1:6)] #choose the level of analysis: currently is ipc_class (subclass and maingroup resulted in excessive combinatorial possibilities)
data2019 <- unique(data2019[,1:6]) #remove duplicated rows (ocurrs when they differ in a ipc level below the one under analysis)
colnames(data2019)[6] <- "ipc_code"

#since this is a knowledge recombination study, remove cases with no recombination (i.e. only one ipc_code):
data2019 <- data2019 %>% 
  group_by(id, pub_date, applicants) %>% 
  filter(n()>1) 

data2019$id <- as.factor(data2019$id)
data2019$ipc_code <- as.factor(data2019$ipc_code)

##Use split to create a list with transaction (i.e. patent) and product categories (i.e. ipc), and turn the result into a transactions object
basket2019 <- split(data2019$ipc_code, data2019$id) %>% 
  as("transactions")

#get maximum no. ipc codes in any patent, to use it as parameter below:
max_ipc_comb <- data2019 %>% 
  group_by(id) %>% 
  summarise(n = n()) %>% 
  summarise(max = max(n))

#Use the transaction object to run association rule mining function parameterised to obtain frequent itemset rather than rules:
#   Notes: min sup is set to lowest value possible so that all itemsets and corresp. frequencies are obtained. Also, length and confidence restrictions are minimised
itemsets2019 <- arules::apriori(basket2019, parameter = list(sup = 0.00001, conf = 0.01, maxlen = max_ipc_comb[[1]], target = "frequent itemsets"))
df_itemsets2019 <- as(itemsets2019, "data.frame")
df_itemsets2019$no_items <- str_count(df_itemsets2019$items, ",") + 1
#===End STEP 2=================================================================================================================================================================================




#===STEP 3 (construct novelty measure) ========================================================================================================================================================
#New dataframe with all ipc codes of a patent in a character string separated by commas:
data2019_collapsed <- data2019 %>% 
  group_by(id) %>% 
  summarise(ipc_codes = paste(ipc_code, collapse = ","))

data2019_collapsed$no_items <- str_count(data2019_collapsed$ipc_codes, ",") + 1 #column with ipc codes counts

#Count all ocurrences of the itemsets of each patent (this is an absolute measure of frequency of the itemset):
data2019_collapsed$abs_frequency <- NA
pb <- txtProgressBar(min = 0, max = nrow(data2019_collapsed), initial = 0, style = 2) #progress bar
for (i in 1:nrow(data2019_collapsed)) { #loop through each patent
  list_of_items <- unlist(strsplit(as.character(data2019_collapsed[i,2]), split=",")) #build a vector with ipc codes
  number_of_items <- as.numeric(data2019_collapsed[i,3])
  focal_count <- invisible(data.frame(inspect(subset(itemsets2019, subset = (items %ain% list_of_items)))))
  focal_count <- sum(focal_count$count, na.rm = T) #sum all ocurrences of the itemset (including ocurrences within larger itemsets)
  data2019_collapsed[i, 4] <- focal_count
  setTxtProgressBar(pb,i) #update progress bar
}
data2019_collapsed <- data.frame(data2019_collapsed)

#To control for the size of the itemset, create a measure of relative frequency of the itemset based on comparisons with itemsets of similar size: (calculating avg abs freq of itemsets of same...
   # sizes and computing ratio of abs freq over the calculated avg):
data2019_collapsed_copy <- data2019_collapsed
data2019_collapsed_copy <- data2019_collapsed_copy %>% group_by(no_items) %>% mutate(group_avg_abs_frequency = mean(abs_frequency))

#The actual measure. It's the individual freq. count in relation to mean freq. count of the patents with same number of items. The lower the more novel:
data2019_collapsed_copy$rel_frequency <- data2019_collapsed_copy$abs_frequency / data2019_collapsed_copy$group_avg_abs_frequency
data2019_collapsed_copy <- data.frame(data2019_collapsed_copy)

#Integrate data needed for development of collaboration variable:
data2019_applicants <- data2019 %>% 
  distinct(id, pub_date, country, applicants, inventors, .keep_all= T)

data2019_collapsed_copy <- merge(data2019_collapsed_copy, data2019_applicants[,-6], by = "id")
#===End STEP 3 =============================================================================================================================================================================




#===STEP 4 (Clean and format applicants and inventors in preparation to collaboration classification and counting) =========================================================================
#loop through each patents, clean up the applicants list by removing the inventors when they are also in applicants, and removing double counting the same applicant, thereby leaving pure... 
      #distinct applicants orgs so I can truly count collaborators (unless applicants = inventors):
data2019_collapsed_copy$applicants_cleaned <- NA
for (j in 1:nrow(data2019_collapsed_copy)) {
  applicants_lst <- strsplit(data2019_collapsed_copy[j,9], "; ")[[1]]
  inventors_lst <- strsplit(data2019_collapsed_copy[j,10], "; ")[[1]]
  applicants_lst <- applicants_lst[!(applicants_lst %in% inventors_lst)] #leave only names that are not in inventors list
  #this is to check if there are two applicants but they're "almost" the same (sometimes is the same but differ in a comma or a character, so this checks if the diff is 1 character, then leave just 1 name)
  if(length(applicants_lst) == 2) {
    a <- tolower(applicants_lst[[1]])
    b <- tolower(applicants_lst[[2]])
    diffe <- mapply(function(x,y) sum(x!=y),strsplit(a,""),strsplit(b,"")) #count different characters
    if(diffe == 1) { applicants_lst <- applicants_lst[[1]] }
  }
  # if there's only one applicant but there's a ";", remove it:
  if(length(applicants_lst) == 1) {
    if(!is.na(applicants_lst)) { #step needed to avoid errors
      if(str_detect(applicants_lst, ";")==T) { applicants_lst <- str_remove(applicants_lst, ";")  }  
    }
  }
  resu <- ifelse(length(applicants_lst)!=0, paste(applicants_lst, collapse = "; "), data2019_collapsed_copy[j,9]) # if all the applicants are inventors, don't leave it in blank, leave applicants as is
  data2019_collapsed_copy[j, 11] <- ifelse(resu=="NA", data2019_collapsed_copy[j,9], resu) #When applicants have an erroneous character, the process returns "NA". Leave the applicants as is if this happens.
}
#===End STEP 4 =================================================================================================================================================================================




#===STEP 5 (Create the collaboration variable) =================================================================================================================================================
data2019_collapsed_copy$no_applicants <- str_count(data2019_collapsed_copy$applicants_cleaned, ";") + 1
#===End STEP 5 ================================================================================================================================================================================




#===STEP 6 (Run ANOVA) =======================================================================================================================================================================
#First, test if the constructed variable is normal (visual test, because Shapiro-Wilk test has a limit of N=5000, and Kolmogorov-Smirnov "becomes highly sensitive" with N>1000):
hist(data2019_collapsed_copy$rel_frequency) #Test one: histogram - seems non-normal
ggqqplot(data2019_collapsed_copy$rel_frequency) #Test two: qqplot - yes, it's non-normal

#Second, exploration:
ggplot(data2019_collapsed_copy, aes(x=no_applicants, y = rel_frequency)) +
  geom_point() +
  theme_bw()

#Because of non normality, I should use non-parametric methods:
cor(data2019_collapsed_copy$no_applicants, data2019_collapsed_copy$rel_frequency, method = "spearman")

#Analysis/method 1: Wilcoxon Signed-Rank Test (two continuous variables):
wilcox.test(data2019_collapsed_copy$no_applicants, data2019_collapsed_copy$rel_frequency, paired = T) #RESULT: p<0.001...collaboration (no_applicants) affects novelty (rel. frequency of IPC class recombinations)

#Analysis/method 2: Mann-Whitney-Wilcoxon Test (one continuous variable and one categorical variable):
data2019_collapsed_copy$no_applicants_groups <- ifelse(data2019_collapsed_copy$no_applicants > 1, 1, 0)
data2019_collapsed_copy$no_applicants_groups <- factor(data2019_collapsed_copy$no_applicants_groups)
wilcox.test(rel_frequency ~ no_applicants_groups, data = data2019_collapsed_copy) #RESULT: p=0.001...collaboration (no_applicants) affects novelty (rel. frequency of IPC class recombinations) WHEN using binary variable for collaboration (0 = no collab, 1 = any collab)

#Analysis/method 3: Kruskal-Wallis Test (one continuous variable and one categorical variable):
data2019_collapsed_copy$no_applicants_groups <- factor(data2019_collapsed_copy$no_applicants)
kruskal.test(rel_frequency ~ no_applicants_groups, data = data2019_collapsed_copy) #RESULT: p<0.0066..collaboration (no_applicants) doesn't affects novelty (rel. frequency of IPC class recombinations)

cutting_collab <- data2019_collapsed_copy %>% 
  group_by(no_applicants) %>% 
  summarise(n = n())

data2019_collapsed_copy$no_applicants_groups <- cut(data2019_collapsed_copy$no_applicants, breaks=c(0,1,2,4,7,19), labels = c("low", "low_mid", "mid", "mid_high", "high")) #form categories cutting in specific values, so that firms with same incorp year are not split into different grous
kruskal.test(rel_frequency ~ no_applicants_groups, data = data2019_collapsed_copy) #RESULT: p<0.00028...collaboration (no_applicants) affects novelty (rel. frequency of IPC class recombinations)

# Perform Dunn Test - post hoc test to detect which specific means are significant from the others
DT <- dunnTest(rel_frequency ~ no_applicants_groups, data = data2019_collapsed_copy, method = "bh") # Adjusts p-values for multiple comparisons;
DT
#Compact letter display:
PT <- DT$res
cldList(P.adj ~ Comparison, data = PT, threshold = 0.05)

#Desc. stats for results table:
desc_stats <- data2019_collapsed_copy %>% group_by(no_applicants_groups) %>% summarise(N = n(), Mean = mean(rel_frequency), Median = median(rel_frequency))
#===End STEP 6 =================================================================================================================================================================================




#===Step 7 (Visualisation) =============================================================================================================================================================
#7A - Top most and least frequent IPC codes:
# Prepare dataset with (individual) item frequency:
df_freq <- df_itemsets2019[df_itemsets2019$no_items==1, 1:2]
df_freq <- df_freq[order(df_freq$support, decreasing = T),]
df_freq$items <- substr(df_freq$items, 2, 4) #remove '{' and '}'
df_hifreq <- head(df_freq, 10)
df_hifreq$group <- "Most frequent items"
df_lofreq <- tail(df_freq, 10)
df_lofreq$group <- "Least frequent items"
df_freq <- rbind(df_hifreq, df_lofreq)
df_freq$items <- factor(df_freq$items, levels = unique(df_freq$items))
df_freq$group <- factor(df_freq$group, levels = unique(df_freq$group))
df_freq$support <- df_freq$support * 100 # to make it a percentage

#Bring in the titles/labels of the IPC codes from db:
ipc_labels <- query("SELECT full_id AS 'items', name FROM CPS_project_patents.ipc_classes;")
df_freq <- merge(df_freq, ipc_labels, by = "items", all = T)
df_freq[df_freq$items=='B68',4] <- "SADDLERY; UPHOLSTERY" #NOTE: I had to add it manually because it's not in the ipc codes tables I've uploaded into mysql db (perhaps they are a few months older than the patent data)
df_freq <- df_freq[complete.cases(df_freq),] #filter out the ipc codes coming from db that weren't in df_freq originally (this is to validate what's missing)
df_freq$name <- ifelse(nchar(df_freq$name) > 50, paste0(substr(df_freq$name, 1, 50), "..."), df_freq$name) 
df_freq$items <- paste(df_freq$items, tolower(df_freq$name), sep = " - ")
df_freq$items <- factor(df_freq$items, levels = df_freq$items[order(df_freq$support)])
#end

ggplot(df_freq, aes(x= items, y = support, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#d95f02", "#1b9e77")) +
  theme_minimal(base_family = "sans", base_size = 22) +
  theme(legend.position = "none", panel.grid = element_blank(), axis.title.x=element_blank(), axis.ticks.y = element_line(), strip.background = element_rect(fill="grey95", colour = "grey90"),
        panel.border = element_rect(colour = "grey90", fill = NA), plot.title = element_text(size=27, face="bold", hjust = 0.5)) +
  coord_flip() +
  facet_wrap(~group, nrow = 2, scales = "free")


#7B - Tree graph (in network representation) with sample of IPC codes combinations:
set.seed(123)

#Prepare dataframe needed for plot (with items and corresponding cuts/levels as columns)
df_itemset_freq <- df_itemsets2019[df_itemsets2019$no_items==4,]
df_itemset_freq <- df_itemset_freq[sample(nrow(df_itemset_freq), 30),] # sample of x random itemsets
df_itemset_freq$items <- gsub('\\{|}', '', df_itemset_freq$items) #remove '{' and '}'
df_itemset_freq$lev1 <- "itemsets"
df_itemset_freq$lev2 <- substr(df_itemset_freq$items,1,7)
df_itemset_freq$lev3 <- substr(df_itemset_freq$items,1,11)
df_itemset_freq$lev4 <- substr(df_itemset_freq$items,1,15)

# create a new column with the all columns collapsed together
columns <- colnames(df_itemset_freq)[5:ncol(df_itemset_freq)]
df_itemset_freq$pathString <- apply(df_itemset_freq[ ,columns] , 1 , paste , collapse = "/" ) 
df_itemset_freq <- df_itemset_freq[,-c(1:4)]

#build the igraph object
g <- data.tree::as.Node(df_itemset_freq)
ig <- data.tree::as.igraph.Node(g)

#build a dataframe (df_2) with key info for the plot such as node's colour and size:
itemset <- unique(c(df_itemset_freq$lev1, df_itemset_freq$lev2, df_itemset_freq$lev3, df_itemset_freq$lev4))
rel_freq <- rep(NA, length(itemset))
df_2 <- data.frame(cbind(itemset, rel_freq))
df_2$itemset <- as.character(df_2$itemset)
df_2$rel_freq <- as.numeric(df_2$rel_freq)
df_itemsets2019_copy <- df_itemsets2019
df_itemsets2019_copy$items <- gsub('\\{|}', '', df_itemsets2019_copy$items) #remove '{' and '}'
pb <- txtProgressBar(min = 0, max = nrow(df_2), initial = 0, style = 1) #progress bar
for(k in 1:nrow(df_2)) { #loop through the unique nodes
  for(l in 1:nrow(df_itemsets2019_copy)) {#loop through the big data table
    if(df_itemsets2019_copy[l,1] == df_2[k,1]) {
      df_2[k,2] <- df_itemsets2019_copy[l,2]
    }
  }
  setTxtProgressBar(pb,k) #update progress bar
}
df_2$rank <- rank(df_2$rel_freq)
nodes_vctr <- names(V(ig))
df_2 <- df_2[match(nodes_vctr, df_2$itemset),] #consistent sorting
V(ig)$color <- df_2$rank
V(ig)$size <- df_2$rank/10
cols <- viridis(100)
pale <- colorRampPalette(cols)
ig$palette <- pale(nrow(df_2))
lay <- layout_with_fr(ig) #type of layout for network plot

#Plot:
par(mar = c(1, 6.5, 1, 6.5) + 1)
plot(ig, layout = lay, vertex.frame.color = "white", vertex.label.color = "black", vertex.label.family="sans", vertex.label.cex = 1.5,
     edge.color = "grey50", edge.width = 0.5, asp = 0.55, margin = -0.17, vertex.label.degree = pi/2, vertex.label.dist = 0)

#Collaboration-novelty effect: (inverted u shape, recalling lower y values means more novelty)
ggplot(data2019_collapsed_copy, aes(x = no_applicants_groups, y = log(rel_frequency), colour = no_applicants_groups)) +
  geom_jitter(alpha = .3, width = .3) +
  geom_boxplot(colour = "grey35", alpha = 0.025) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 1), plot.title = element_text(hjust = 0.5))
#===End STEP 7 =================================================================================================================================================================================

lapply(dbListConnections(dbDriver(drv = "MySQL")), dbDisconnect) # close all opened connections to db



