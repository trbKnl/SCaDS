################## ADNI data  #################################


require(mice)
require(psych)

MyData <- data.frame(read.csv(file="./ADNI/SelectedDataCombined.csv", 
                              header=F, sep=",", stringsAsFactors=FALSE)[, -746])   #note 746th columns contains nothing. the original data have 645 columns, but somehow it has an extra column once load to R.

Data_ADNIgo <- MyData[, MyData[1, ]=="ADNIGO"]

ADNIgo <- data.frame(Data_ADNIgo, stringsAsFactors=FALSE)


RosterID <- array()
for (i in 1:length(SubjectID)){
  RosterID[i] <- strsplit(as.character(SubjectID[i]), split = "_")[[1]][3]
}

ADNIgo_final <- rbind(ADNIgo, RosterID)
write.csv(ADNIgo_final, file = "./ADNI/ADNIgo_final.csv")



################### processing ECOGPT.csv data ####################
ECOGPT <- data.frame(read.csv(file="./ADNI/DataUsedforPaper/ECOGPT.csv", 
                              header=T, sep=",", stringsAsFactors=FALSE))

library(mice)

# Memory subscale
a <- match("MEMORY1", names(ECOGPT))
b <- match("MEMORY8", names(ECOGPT))
MEMORY <- ECOGPT[, a:b]
MEMORY[which(MEMORY==9, arr.ind = T)] <- NA
MEMORY_IMP <- mice(MEMORY, seed = 1,  printFlag=F)
MEMORY_IMP$imp
MEMORY_Final <- complete(MEMORY_IMP)
Sum_Memory <- rowSums(MEMORY_Final)

# Language subscale
a <- match("LANG1", names(ECOGPT))
b <- match("LANG9", names(ECOGPT))
LANG <- ECOGPT[, a:b]
LANG[which(LANG==9, arr.ind = T)] <- NA
LANG_IMP <- mice(LANG, seed = 1,  printFlag=F)
LANG_Final <- complete(LANG_IMP)
Sum_Lang <- rowSums(LANG_Final)

# VISUAL-SPATIAL AND PERCEPTUAL ABILITIES scale
a <- match("VISSPAT1", names(ECOGPT))
b <- match("VISSPAT8", names(ECOGPT))
VISSPAT <- ECOGPT[, a:b]
VISSPAT[which(VISSPAT==9, arr.ind = T)] <- NA
VISSPAT_IMP <- mice(VISSPAT, seed = 1,  printFlag=F)
VISSPAT_Final <- complete(VISSPAT_IMP)
Sum_Visspat <- rowSums(VISSPAT_Final)

# EXECUTIVE FUNCTIONING: PLANNING scale
a <- match("PLAN1", names(ECOGPT))
b <- match("PLAN5", names(ECOGPT))
PLAN <- ECOGPT[, a:b]
PLAN[which(PLAN==9, arr.ind = T)] <- NA
PLAN_IMP <- mice(PLAN, seed = 1,  printFlag=F)
PLAN_Final <- complete(PLAN_IMP)
Sum_Plan <- rowSums(PLAN_Final)

# EXECUTIVE FUNCTIONING: ORGANIZATION scale
a <- match("ORGAN1", names(ECOGPT))
b <- match("ORGAN6", names(ECOGPT))
ORGAN <- ECOGPT[, a:b]
ORGAN[which(ORGAN==9, arr.ind = T)] <- NA
ORGAN_IMP <- mice(ORGAN, seed = 1,  printFlag=F)
ORGAN_Final <- complete(ORGAN_IMP)
Sum_Organ <- rowSums(ORGAN_Final)


# EXECUTIVE FUNCTIONING: DIVIDED ATTENTION scale
a <- match("DIVATT1", names(ECOGPT))
b <- match("DIVATT4", names(ECOGPT))
DIVATT <- ECOGPT[, a:b]
DIVATT[which(DIVATT==9, arr.ind = T)] <- NA
DIVATT_IMP <- mice(DIVATT, seed = 1,  printFlag=F)
DIVATT_Final <- complete(DIVATT_IMP)
Sum_Divatt <- rowSums(DIVATT_Final)

ECOGPT_IMP <- cbind(ECOGPT[, 1:9], MEMORY_Final, LANG_Final, VISSPAT_Final, PLAN_Final, ORGAN_Final, DIVATT_Final, Sum_Memory, Sum_Lang, Sum_Visspat, Sum_Plan, Sum_Organ, Sum_Divatt)
write.csv(ECOGPT_IMP, file = "./ADNI/DataUsedforPaper/ECOGPT_imp.csv")



############################# process the merged data Final_merged.csv ############################

Merged_D <- data.frame(read.csv(file="./ADNI/DataUsedforPaper/merge data/Final_merged.csv", 
                    header=F, sep=",", stringsAsFactors=FALSE))

Merged_D <- data.frame(read.csv(file="./ADNI/DataUsedforPaper/merge data/Final_merged.csv", 
                                header=F, sep=",", stringsAsFactors=FALSE))  #laptop office.

# randomly pick 20 persons 
#set.seed(112)
#index <- sample(2:296, 20)
#Merged_selected <- Merged_D[c(1,index), ]
#Merged_selected <- Merged_selected[c(-3, -4, -5, -9, -14, -16, -17, -20), ]  #Note that here I removed the subjects that contain missing values. Here, the missing values are because they are not measured. 
                                                         #Of course we can do imputation, but I do not think it makes sense here. 
neuropsy <- Merged_D[, 2:13] #the first column is left out, since its subject ID
genes <- Merged_D[, 15:1140] #the first column is left out, since its subject ID

neuropsy_names <- neuropsy[1, ]
neuropsy_data <- data.matrix(neuropsy[2:296, ])
genes_names <- genes[1, ] 
genes_data <- data.matrix(genes[2:296, ])

library(psych)
describe(neuropsy_data) 
# we remove the NA entries in columns V2 - V7
# we remove the -1 entries (missing value) in V8 - V13
table(rowSums(is.na(neuropsy_data))) # 117 enties contain NA's for the entire 12 questions, remove them.
which(rowSums(is.na(neuropsy_data)) != 12)
neuropsy_data <- neuropsy_data[which(rowSums(is.na(neuropsy_data)) != 12), ]
genes_data <- genes_data[which(rowSums(is.na(neuropsy_data)) != 12),]

minus_index <- neuropsy_data[, 12] != -1 #a few rows contain missing values coded as "-1"
neuropsy_data <- neuropsy_data[minus_index, ]  
genes_data <- genes_data [minus_index, ]


# some genes are measured repeatedly. Here we keep the first measure. 
index_toremove <- grepl("_", genes_names)
genes_names <- genes_names[!index_toremove]
genes_data <- genes_data[, !index_toremove]

colnames(neuropsy_data) <-  neuropsy_names
colnames(genes_data) <- genes_names
 


save(neuropsy_data, genes_data, file = "./ADNI/DataUsedforPaper/merge data/ADNI_final.RData")

