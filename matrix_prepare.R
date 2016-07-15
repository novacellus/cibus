#Prepare corpus matrix
library(MAAS)
library(wordspace)
#1. Read-in 
#data <- "~/Kod/RWorkspace/colloc/data/patrologia_bb.entire_lemmas.tbl"
data <- "~/Kod/RWorkspace/wordspace/tables/patrologia.full.tagg.tbl" #Tagging changed
lemma_triples <- as.data.table(read.csv(data,
                          sep = "\t",header = F))
colnames(lemma_triples) <- c("freq","word1", "word2")

#2. Constructing matrix
lemma_VObj <- dsm(target = lemma_triples$word1,
                       feature = lemma_triples$word2,
                       score = lemma_triples$freq,
                       raw.freq=TRUE, sort=TRUE,verbose = T)
dim(lemma_VObj)

#3. Reducing matrix
#a. Cutting out low-frequency words
lemma_VObj <- subset(lemma_VObj, nnzero >= 3, nnzero >= 3, recursive=TRUE)

#b. Weighthing of Matrix
#Parameters:
#@score: simple-ll, ?
#@transform: log, ?
#@normalize: TRUE
#@scale: standardize (z-transformation), center, scale
lemma_VObj_weight <- dsm.score(lemma_VObj, score="simple-ll", transform="log", normalize=TRUE)

#c. Reduce by subspace projection
#@method: svd, rsvdm, asvd, ri, ri+rsvd
#@n: target dimensions
#@oversampling: default - 2 for rsvd, 10 for other (>7 mins for Patrologia)
lemma_VObj_proj  <- dsm.projection(lemma_VObj_weight, method="rsvd", n=300, oversampling=4)
save(lemma_VObj_proj,file = "lemma_VObj_proj")

load("lemma_VObj_proj")
