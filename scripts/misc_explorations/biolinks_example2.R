library("TCGAbiolinks")
query <- GDCquery(   project = "TCGA-LIHC",  
                     data.category = "Clinical", 
                     data.type = "Clinical Supplement",
                     dat.format = "BCR Biotab",
                     legacy = FALSE)
GDCdownload(query,directory = ".")
clinical <- GDCprepare(query,directory = ".")






query <- GDCquery(project = "TCGA-LIHC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
followup <- clinical.BCRtab.all$clinical_follow_up_v4.0_lihc
