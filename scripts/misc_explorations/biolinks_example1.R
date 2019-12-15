library("TCGAbiolinks")
library("tidyverse")
library("magrittr")


#getting xena
xena <- read_tsv("xena_lihc_surv.tsv")
#xml query
query <- GDCquery(   project = "TCGA-LIHC",  
                     data.category = "Clinical", 
                     file.type = "xml", 
                     legacy = FALSE)
GDCdownload(query,directory = ".")
#getting xml
clinical <- GDCprepare_clinic(query, clinical.info = "patient",directory = ".")
survival_data <- as_tibble(clinical[,c("days_to_last_followup","days_to_death","vital_status","bcr_patient_barcode")]) 
survival_data <- filter(survival_data,!is.na(days_to_last_followup)|!is.na(days_to_death))  #not both NA
survival_data$days_to_last_followup[is.na(survival_data$days_to_last_followup)] <- -Inf
survival_data$days_to_death[is.na(survival_data$days_to_death)] <- -Inf
survival_data$os.time<-  pmax(survival_data$days_to_last_followup,survival_data$days_to_death)
survival_data <- filter(survival_data,os.time>0 & !is.na(os.time)) #ensuring positive values
dim(survival_data) #there should be 371 rows here! but there are 343!

clinical2 <- GDCprepare_clinic(query, clinical.info = "follow_up",directory = ".")
survival_data2 <- as_tibble(clinical2[,c("days_to_last_followup","days_to_death","vital_status","bcr_patient_barcode")]) 
survival_data2 <- filter(survival_data2,!is.na(days_to_last_followup)|!is.na(days_to_death))  #not both NA
survival_data2$days_to_last_followup[is.na(survival_data2$days_to_last_followup)] <- -Inf
survival_data2$days_to_death[is.na(survival_data2$days_to_death)] <- -Inf
survival_data2$os.time<-  pmax(survival_data2$days_to_last_followup,survival_data2$days_to_death)
survival_data2 <- filter(survival_data2,os.time>0& !is.na(os.time) ) #ensuring positive values
bound <- rbind(survival_data,survival_data2)
bound <- arrange(bound,bcr_patient_barcode,desc(os.time))
bound <- bound[!duplicated(bound$bcr_patient_barcode),] 
bound <- bound[,colnames(xena)]

#getting indexed (directly from gdcquery_clinic)
index.clinical <- GDCquery_clinic("TCGA-LIHC",type = "clinical")
index.clinical <- index.clinical[,c("days_to_last_follow_up","days_to_death","vital_status","submitter_id")]
colnames(index.clinical) <-  c("days_to_last_followup","days_to_death","vital_status","bcr_patient_barcode")
index.clinical <- filter(index.clinical,!is.na(days_to_last_followup)|!is.na(days_to_death))  #not both NA
index.clinical$days_to_last_followup[is.na(index.clinical$days_to_last_followup)] <- -Inf
index.clinical$days_to_death[is.na(index.clinical$days_to_death)] <- -Inf
index.clinical$os.time <- pmax(index.clinical$days_to_last_followup,index.clinical$days_to_death)
index.clinical <- filter(index.clinical,os.time>0 & !is.na(os.time) ) #ensuring positive values
index.clinical <- index.clinical[!duplicated(index.clinical$bcr_patient_barcode),] 
index.clinical <- arrange(index.clinical,bcr_patient_barcode,desc(os.time))
index.clinical <- bound[!duplicated(index.clinical$bcr_patient_barcode),] 
index.clinical <- index.clinical[,colnames(xena)]


#getting biotab
#biotab_query
query <- GDCquery(project = "TCGA-LIHC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
followup_biotab <- clinical.BCRtab.all$clinical_patient_lihc[c("last_contact_days_to","death_days_to","vital_status","bcr_patient_barcode")]
colnames(followup_biotab) <- c("days_to_last_followup","days_to_death","vital_status","bcr_patient_barcode")
followup_biotab <- followup_biotab[3:dim(followup_biotab)[1],]
followup_biotab %<>% mutate_all(~ replace(., . == "[Not Available]" | . == "[Not Applicable]" , NA))
followup_biotab$days_to_death %<>% as.numeric()
followup_biotab$days_to_last_followup %<>% as.numeric()
followup_biotab <- filter(followup_biotab,!is.na(days_to_last_followup)|!is.na(days_to_death))  #not both NA
followup_biotab$days_to_last_followup[is.na(followup_biotab$days_to_last_followup)] <- -Inf
followup_biotab$days_to_death[is.na(followup_biotab$days_to_death)] <- -Inf
followup_biotab$os.time <- pmax(followup_biotab$days_to_last_followup,followup_biotab$days_to_death)
followup_biotab <- filter(followup_biotab,os.time>0 & !is.na(os.time)) #ensuring positive values
followup_biotab2 <- clinical.BCRtab.all$clinical_follow_up_v4.0_lihc[,c("last_contact_days_to","death_days_to","vital_status","bcr_patient_barcode")]
colnames(followup_biotab2) <- c("days_to_last_followup","days_to_death","vital_status","bcr_patient_barcode")
followup_biotab2 <- followup_biotab2[3:dim(followup_biotab2)[1],]
followup_biotab2 %<>% mutate_all(~ replace(., . == "[Not Available]" | . == "[Not Applicable]" , NA))
followup_biotab2$days_to_death %<>% as.numeric()
followup_biotab2$days_to_last_followup %<>% as.numeric()
followup_biotab2 <- filter(followup_biotab2,!is.na(days_to_last_followup)|!is.na(days_to_death))  #not both NA
followup_biotab2$days_to_last_followup[is.na(followup_biotab2$days_to_last_followup)] <- -Inf
followup_biotab2$days_to_death[is.na(followup_biotab2$days_to_death)] <- -Inf
followup_biotab2$os.time <- pmax(followup_biotab2$days_to_last_followup,followup_biotab2$days_to_death)
followup_biotab2 <- filter(followup_biotab2,os.time>0 & !is.na(os.time)) #ensuring positive values
bound2 <- rbind(followup_biotab,followup_biotab2)
bound2 <- arrange(bound2,bcr_patient_barcode,desc(os.time))
bound2 <- bound2[!duplicated(bound2$bcr_patient_barcode),] 
bound2 <- bound2[,colnames(xena)]



#checking_equality
all(index.clinical == bound)
all(bound2 == bound)
all(index.clinical == bound2)
all(bound == xena) #have some inequalities

dif_ostime <- bound$os.time != xena$os.time
dif_barcode <- bound$bcr_patient_barcode != xena$bcr_patient_barcode
dif_vital <- bound$vital_status != xena$vital_status

xena_bound_bound <- cbind(xena,bound)
xena_bound_bound <- xena_bound_bound[(dif_ostime | dif_barcode | dif_vital),]
print(xena_bound_bound)
