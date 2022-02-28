#Le librerie necessarie per realizzare gli esercizi sono:
#Comandi base di R
#per più approfondimenti: R for Beginners, Emmanuel Paradis, https://cran.r-project.org/doc/contrib/Paradis-rdebuts_en.pdf




#Lezione 3: 
#(Ref:Finding Candidate Binding Sites for Known Transcription Factors via Sequence Matching https://www.bioconductor.org/packages/release/workflows/vignettes/generegulation/inst/doc/generegulation.html)

library(BiocManager)
BiocManager::install(c("MotifDb",  "GenomicFeatures", 
                       "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene",
                       "org.Sc.sgd.db", "BSgenome.Scerevisiae.UCSC.sacCer3",
                       "motifStack", "seqLogo"))
  
  