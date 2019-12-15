query <- GDCquery(project = "TCGA-CESC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", #this is necessary
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)

TCGAbiolinks:::getProjectSummary("TCGA-CESC")
getGDCprojects()$project_id


#used indexed data for followup 
# use the clinical supplement for hpv


patient_cesc <- clinical.BCRtab.all$clinical_patient_cesc