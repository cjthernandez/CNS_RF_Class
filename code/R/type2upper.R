### Convert tumor classes to tumor families (As defined by Capper et al. (2018) Figure 1)

# From methylation class get methylation family
type2upper <- function(samp_vec, p=F, met=F){
  # Methylation class families
  other <- c("CHGL", "LGG, SEGA", "LGG, PA PF", "LGG, PA MID", "ANA PA", "HGNET, MN1", "IHG", "LGG, MYB",
             "LGG, PA/GG ST", "PXA")
  nerve <- c("SCHW", "SCHW, MEL")
  pineal <- c("PTPR, A", "PTPR, B", "PIN T,  PB A", "PIN T,  PB B", "PIN T, PPT")
  mesenc <- c("CHORDM", "EWS", "HMB", "MNG", "SFT HMPC", "EFT, CIC")
  if(!met){
    melano <- c("MELCYT", "MELAN")
  }else{
    melano <- c("MELCYT")
  }
  plexus <- c("PLEX, AD", "PLEX, PED A", "PLEX, PED B")
  glioIDH <- c("A IDH", "A IDH, HG", "O IDH")
  hemat <- c("LYMPHO", "PLASMA")
  epend <- c("EPN, RELA", "EPN, PF A", "EPN, PF B", "EPN, YAP", "EPN, SPINE", "EPN, MPE", "SUBEPN, PF",
             "SUBEPN, ST", "SUBEPN, SPINE")
  sella <- c("CPH, ADM", "CPH, PAP", "PITAD, ACTH", "PITAD, STH DNS B", "PITAD, PRL", "PITAD, FSH LH",
             "PITAD, STH SPA", "PITAD, STH DNS A", "PITAD, TSH", "PITUI")
  embryo <- c("ETMR", "MB, WNT", "MB, G3", "MB, G4", "MB, SHH INF", "MB, SHH CHL AD", "ATRT, SHH",
              "ATRT, TYR", "ATRT, MYC", "CNS NB, FOXR2", "HGNET, BCOR")
  gliob <- c("DMG, K27", "GBM, G34", "GBM, MES", "GBM, RTK II", "GBM, RTK III", "GBM, RTK I", "GBM, MYCN", "GBM, MID")
  glion <- c("CN", "DLGNT", "LIPN", "LGG, DIG/DIA", "LGG, DNT", "LGG, RGNT", "LGG, GG", "RETB", "ENB, A", "ENB, B", "PGG, nC")
  control <- c("CONTR, ADENOPIT", "CONTR, WM", "CONTR, CEBM", "CONTR, HEMI", "CONTR, HYPTHAL", "CONTR, INFLAM",
               "CONTR, PINEAL", "CONTR, PONS", "CONTR, REACT")
  # Classes from GSE44661 and GSE108576
  metastasis <- c("Breast cancer brain metastasis","Lung cancer brain metastasis","Melanoma brain metastasis",
                  "Uncertain primary tumor brain metastasis")
  
  classes <- c(control, other, nerve, pineal, mesenc, melano, plexus, glioIDH, hemat, epend, sella, glion, gliob, embryo, metastasis)
  classes_point <- gsub("[^[:alnum:]]",".",classes)
  
  class_fam <- c(rep("Control", length(control)),
                 rep("Other glioma", length(other)),
                 rep("Nerve", length(nerve)),
                 rep("Pineal", length(pineal)),
                 rep("Mesenchymal", length(mesenc)),
                 rep("Melanocytic", length(melano)),
                 rep("Plexus", length(plexus)),
                 rep("Glioma IDH", length(glioIDH)),
                 rep("Hematopoietic", length(hemat)),
                 rep("Ependymal", length(epend)),
                 rep("Sella", length(sella)),
                 rep("Glio-neuronal", length(glion)),
                 rep("Glioblastoma", length(gliob)),
                 rep("Embryonal", length(embryo)),
                 rep("Brain metastasis", length(metastasis)))
  
  up_clust <- rep("0", length(samp_vec))
  
  ii <- match(samp_vec, classes)
  
  jj <- which(!is.na(ii))

  up_clust[jj] <- class_fam[na.omit(ii)]
  
  up_clust[grep("New", samp_vec)] <- "New sample"
  
  return(up_clust)
}

# From methylation family get list of methylation classes
upper2type <- function(upper, met=F){
  other <- c("CHGL", "LGG, SEGA", "LGG, PA PF", "LGG, PA MID", "ANA PA", "HGNET, MN1", "IHG", "LGG, MYB",
             "LGG, PA/GG ST", "PXA")
  nerve <- c("SCHW", "SCHW, MEL")
  pineal <- c("PTPR, A", "PTPR, B", "PIN T,  PB A", "PIN T,  PB B", "PIN T, PPT")
  mesenc <- c("CHORDM", "EWS", "HMB", "MNG", "SFT HMPC", "EFT, CIC")
  if(!met){
    melano <- c("MELCYT", "MELAN")
  }else{
    melano <- c("MELCYT")
  }
  plexus <- c("PLEX, AD", "PLEX, PED A", "PLEX, PED B")
  glioIDH <- c("A IDH", "A IDH, HG", "O IDH")
  hemat <- c("LYMPHO", "PLASMA")
  epend <- c("EPN, RELA", "EPN, PF A", "EPN, PF B", "EPN, YAP", "EPN, SPINE", "EPN, MPE", "SUBEPN, PF",
             "SUBEPN, ST", "SUBEPN, SPINE")
  sella <- c("CPH, ADM", "CPH, PAP", "PITAD, ACTH", "PITAD, STH DNS B", "PITAD, PRL", "PITAD, FSH LH",
             "PITAD, STH SPA", "PITAD, STH DNS A", "PITAD, TSH", "PITUI")
  embryo <- c("ETMR", "MB, WNT", "MB, G3", "MB, G4", "MB, SHH INF", "MB, SHH CHL AD", "ATRT, SHH",
              "ATRT, TYR", "ATRT, MYC", "CNS NB, FOXR2", "HGNET, BCOR")
  gliob <- c("DMG, K27", "GBM, G34", "GBM, MES", "GBM, RTK II", "GBM, RTK III", "GBM, RTK I", "GBM, MYCN", "GBM, MID")
  glion <- c("CN", "DLGNT", "LIPN", "LGG, DIG/DIA", "LGG, DNT", "LGG, RGNT", "LGG, GG", "RETB", "ENB, A", "ENB, B", "PGG, nC")
  control <- c("CONTR, ADENOPIT", "CONTR, WM", "CONTR, CEBM", "CONTR, HEMI", "CONTR, HYPTHAL", "CONTR, INFLAM",
               "CONTR, PINEAL", "CONTR, PONS", "CONTR, REACT")
  # Classes from GSE44661 and GSE108576
  metastasis <- c("Breast cancer brain metastasis","Lung cancer brain metastasis","Melanoma brain metastasis")#,
                  #"Uncertain primary tumor brain metastasis") ## UNCERTAIN REMOVED
  
  switch(upper,
         "Control"=control,
         "Other glioma"=other,
         "Nerve"=nerve,
         "Pineal"=pineal,
         "Mesenchymal"=mesenc,
         "Melanocytic"=melano,
         "Plexus"=plexus,
         "Glioma IDH"=glioIDH,
         "Hematopoietic"=hemat, 
         "Ependymal"=epend, 
         "Sella"=sella, 
         "Glio-neuronal"=glion,
         "Glioblastoma"=gliob,
         "Embryonal"=embryo,
         "Brain metastasis"=metastasis)
}