library("TCGAbiolinks")
library("tidyverse")
library("magrittr")



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



