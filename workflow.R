#Dictionary similarity

library(tm)
library(lsa)
library(koRpus)
library(reshape2)
library(data.table)
library(magrittr)
#Wczytaj slownik
lewis <- as.data.frame(read.csv(file="dicts/Lewis2table.csv"))
defs <- stringr::str_replace_all(lewis$defs,pattern = "[˘¯]","")
write.table(defs,file="defs.txt",quote = FALSE,row.names = FALSE,col.names = F)
#Usunac wiersz tytulowy!!!
tagged <- treetag("defs.txt", treetagger = "manual",lang = "en", 
                  TT.options  = list(path="/home/krzys/Programy/treetagger", 
                                     tokenizer="utf8-tokenize.perl", 
                                     tagger="tree-tagger",
                                     params="english-utf8.par",
                                     abbrev="english-abbreviations",
                                     no.unknown=TRUE),
                  stopwords=tm::stopwords("en"),stemmer = SnowballC::wordStem)
#Oblicza liczbe wyrazow w wierszu
counts <- as.data.frame(matrix(nrow=length(stringr::str_count(defs,stringr::boundary("word"))),ncol = 4))
colnames(counts) <- c("prev","pres","start","end")

i = 0
prev_count = 0
for (pres_count in stringr::str_count(defs,stringr::boundary("word") ) ) {
  i = i+1
  if (pres_count == 0) {
    prev_count <- prev_count + pres_count
    a <- list(pres_count,prev_count,0,0)
  } else {
    start = prev_count + 1
    end = start + pres_count - 1
    a <- list(prev_count,pres_count,start,end)
  }
  counts[i,] <- a
  prev_count <- prev_count + pres_count
}
rm(list=c("i","a","start","end","pres_count","prev_count"))

lemmas_w_tags <- paste(tagged@TT.res$lemma,tagged@TT.res$tag, sep = "_")
lemmas_w_tags_w_stop_tmp <- paste (lemmas_w_tags, 
                                   as.character(tagged@TT.res$stop), 
                                   sep = "_")
lemmas_w_tags <- apply(X = counts,1,FUN = function(row) paste ( lemmas_w_tags[row[["start"]]:row[["end"]] ], collapse = " ") )
lemmas_w_tags_w_stop<- apply(X = counts,1,FUN = function(row) paste ( lemmas_w_tags_w_stop_tmp[row[["start"]]:row[["end"]] ], collapse = " ") )

lemmas_w_tags_nostop <- apply(X = counts,1,FUN = function(row) paste ( lemmas_w_tags_w_stop_tmp[row[["start"]]:row[["end"]] ][ grep(x = lemmas_w_tags_w_stop_tmp[row[["start"]]:row[["end"]] ], pattern = "_FALSE") ], collapse = " ") )
lemmas_w_tags_nostop_wo_labels <- stringr::str_replace_all(lemmas_w_tags_nostop, pattern = "_FALSE",replacement = "")
#lemmas_w_tags_nostop <- apply(X = counts,1,FUN = function(row) paste ( lemmas_w_tags_w_stop[row[["start"]]:row[["end"]]][grep(x = lemmas_w_tags_w_stop[row[["start"]]:row[["end"]] ], pattern = "_FALSE")], collapse = " ") )


lemmas <- apply(X = counts,1,FUN = function(row) paste ( tagged@TT.res$lemma[row[["start"]]:row[["end"]] ], collapse = " ") )
wclasses <- apply(X = counts,1,FUN = function(row) paste ( tagged@TT.res$wclass[row[["start"]]:row[["end"]] ], collapse = " ") )
tags <- apply(X = counts,1,FUN = function(row) paste ( tagged@TT.res$tag[row[["start"]]:row[["end"]] ], collapse = " ") )

annotation <- data.frame(def_lemmas=lemmas,def_tags=tags,def_wclasses=wclasses,def_lemmas_w_tags=lemmas_w_tags,no_stops = lemmas_w_tags_nostop_wo_labels)

lewis_annotated <- cbind(lewis,annotation)

rm(list=c("counts", "defs","tagged",
          "lemmas", "lemmas_w_tags",
          "lemmas_w_tags_w_stop","lemmas_w_tags_no_stop",
          "lemmas_w_tags_nostop", "lemmas_w_tags_nostop_wo_labels",
          "lemmas_w_tags_w_stop_tmp",
          "lemmas_w_tags_no_stop_wo_labels",
          "wclasses","tags","annotation"))

#Extract list of glosses
no_stops_levels <- lewis_annotated$no_stops  %>% 
  sapply(.,stringr::str_split,pattern=" ") %>% 
  unlist %>% unique
gloss_nostop <-list()
gloss_nostop <- apply(lewis_annotated,1,
                      function(row) gloss_nostop[[as.character(row["lemma"])]] <- 
                        unlist(as.vector(stringr::str_split(row["no_stops"],pattern = " "))) )
names(gloss_nostop) <- lewis_annotated$lemma

rm("no_stops_levels")

gloss_nostop_melted <- data.table(melt(gloss_nostop))
rm("gloss_nostop")

gloss_nostop_wide <-dcast(gloss_nostop_melted,
                           formula= L1  ~ value )

#Poszukiwanie wyrazu w definicji
gloss_nostop_wide[which(gloss_nostop_wide[,food_NN > 0])]$L1




#Obliczam iloczyn każdego wektora z każdym
#gloss_inter <- lapply(gloss_nostop, function(vector_this) {
#  lapply(gloss_nostop, function(vector_that) {
#    length ( intersect(vector_this, vector_that) )
#  })
#})


#Test

#gloss_inter <- lapply(gloss_nostop[1:100], function(vector_this) {
#  lapply(gloss_nostop, function(vector_that) {
#    length ( intersect(vector_this, vector_that) )
#  })
#})

#gloss_df <- as.data.frame(matrix(nrow=length(gloss_nostop),ncol=length(gloss_nostop)))
#colnames(gloss_matrix) <- names(gloss_nostop)
#rownames(gloss_matrix) <- names(gloss_nostop)

#i = 0
#gloss_inter_mc <- lapply(gloss_nostop[2], function(vector_this) {
#  lapply(gloss_nostop, function(vector_that) {
#    i <- i+1
#    gloss_matrix[i,] <- length ( intersect(vector_this, vector_that) )
#  })
#})

#lapply(gloss_nostop[2], function(vector_this) {
#  lapply(gloss_nostop, function(vector_that) {
#  gloss_matrix[1,] <- length ( intersect(vector_this, vector_that) )
#  })
#})




#DocumentTermMatrix
#corp <- tm::VCorpus(x = VectorSource(lewis_annotated$no_stops) )
#meta(corp,type = "local",tag = "lemma") <- lewis_annotated$lemma
#meta(corp,type = "local",tag = "id") <- lewis_annotated$lemma

#corp_dtm <- DocumentTermMatrix(corp)
#corp_tdm <- TermDocumentMatrix(corp)

#Terminy wystepujace min. 150 razy
#findFreqTerms(corp_dtm,150)

#dtm_lsa <- as.data.frame(inspect(corp_tdm))
#rownames(dtm_lsa) <- corp_dtm$dimnames$Terms
#colnames(dtm_lsa) <- lewis_annotated$lemma
