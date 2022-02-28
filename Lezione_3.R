#Predizione in silico delle regioni consenso per fattori di trascrizione conosciuti (C/EBP??)
#del gene TP53 
#seguendo l'esempio riportato su 
#https://www.bioconductor.org/packages/release/workflows/vignettes/generegulation/inst/doc/generegulation.html

#STEP1: Library activation
library(MotifDb)
library(S4Vectors)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Sc.sgd.db)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

#STEP2: Creare un *.txt (in questo caso chiamato TP53sequenceFASTA.txt) con la 
#WHOLE SEQUENCE del gene TP53 in fasta (https://www.ncbi.nlm.nih.gov/nuccore/NG_017013.2?from=5001&to=24149&report=fasta)
#e quindi nomina il tuo gene come un type character [1:n] usando la scansione dei dati sul *.txt  
#ATTENZIONE: trasformare il testo nel *.txt in una sola riga con la sequenza da analizzare. Per farlo:  
TP53raw<- data.frame(scan("TP53sequenceFASTA.txt", what="character")) #importa il *.txt come dataframe
TP53singleline <- c(paste(TP53raw$scan..TP53sequenceFASTA.txt...what....character.., collapse= "")) #converto il dataframe in una sola riga con la sequenza di TP53
TP53

#STEP2.1: cercare se il fattore di trascrizione di nostro interesse è presente nel 
#pacchetto "MotifDb"

MotifDb [grep ('CEBPB', values (MotifDb)$geneSymbol, ignore.case=TRUE)]
Hsapiens_jaspar2018_CEBPB <- query(MotifDb, "CEBPB") [[1]]  
Hsapiens_jaspar2018_CEBPB_nouno <- query(MotifDb, "CEBPB")

#N.B.: se non sai qual'è la nomenclatura (official symbol) del tuo gene prova a cercarlo 
#sul genebase del NCBI. Nel nostro caso la nomenclatura di C/EBP?? è CEBPB.

#N.B.1: su wikipedia puoi trovare una lista di tutti i fattori di trascrizione umani 
#conosciuti https://en.wikipedia.org/wiki/List_of_human_transcription_factors

#N.B.2: se non trovi il fattore di trascrizione, assicurati che nel database usato non abbia una
#nomenclatura vecchia. Ad esempio il fattore di trascrizione GTF3A, è chiamato AP2 nel nostro 
#database, che è il nome del gene omologo nella specie vegetale chiamata Arabidopsis thaliana.
#Infatti se cerchiamo:
#MotifDb [grep ('GTF3A', values (MotifDb)$geneSymbol, ignore.case=TRUE)]
#l'oggetto risultante avrà valore 0. Ma se cerchiamo:
#MotifDb [grep ('AP2', values (MotifDb)$geneSymbol, ignore.case=TRUE)]
#avremo un motifDb di lunghezza 122, con agiornamento al 2013:
#| Created from downloaded public sources: 2013-Aug-30
#| 122 position frequency matrices from 14 sources:

#STEP3: trovare qual'è il motif necessario per trovare il sito di legame tra la regione consenso e CEBPB
seqLogo(Hsapiens_jaspar2018_CEBPB)

#appare un plot con il motif più probabile, in questo caso ATTACGCAA. Pattern che corrisponde a
#quello riportato in #bibliografia RTTGCGYAAY (R = A or G, and Y = C or T) come riportato 
#nell'articolo https://www.jbc.org/article/S0021-9258(17)45942-4/fulltext

#STEP4: matchare la sequenza TP53 con il motif Hsapiens_jaspar2018_CEBPB per cercare il sito 
#di legame del fattore di trascrizione di riferimento (CEBPB) 
pcm.CEBPB <- round(100 * Hsapiens_jaspar2018_CEBPB)

#match del promoter con la nostra sequenza di DNA
plus_hits <- matchPWM(pcm.CEBPB, TP53singleline [[1]], "90%") #La percentuale di match dovrebbe essere maggiore del 90%
plus_hits
#STEP5: crea un dataframe per poter salvare il risultato in un file txt

exporting <- data.frame(TP53singleline [[1]],
                              start=start(plus_hits), end=end(plus_hits),
                              strand="+",
                              seq=as.character(plus_hits))
#togli la prima colonna, o sul risultato ti ritroverai tutta la sequenza di TP53 
exporting$TP53singleline..1.. <- NULL

write.csv(exporting, file = "Results/Lezione3_Result.txt")

#per salvare il plot del motif
dev.copy(tiff,'Results/Lezione_3_plot_motifregioneconsenso.tiff')
dev.off()
