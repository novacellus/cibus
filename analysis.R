library(igraph)
library(wordspace)
library(dplyr)
##Dictionary
#Extract words from LS dictionary
food_LS <- gloss_nostop_wide[which(gloss_nostop_wide[,food_NN > 0])]$L1
#Annotated with PoS
food_LS <- read.csv("export/food_LS_annot.csv",
                    header = T,sep = ",")[,2:3]
#Generate matrix of definitions for selected words only
food_LS_mat <- gloss_nostop_wide[food_LS]
food_LS_MAT_null <- 
  apply(food_LS_mat, 2, function(x) all( x == 0 )) #Search form empty columns and save their names
food_LS_mat <- food_LS_mat[, .SD, 
                           .SDcols=which(food_LS_MAT_null == FALSE) ]
#Clustering on definition matrix
food_LS_mat.df <- as.data.frame(food_LS_mat) #Generates df (only for labelling)
food_LS_mat.df$pos <- NULL
food_LS_mat.df$L1 <- NULL
rownames(food_LS_mat.df) <- food_LS_mat$L1
#The same, but binary
food_LS_mat.df.bin <- food_LS_mat.df
food_LS_mat.df.bin[food_LS_mat.df > 0] <- 1 #Binary

#Long table for visual in gephi
food_LS_mat_long <- food_LS_mat.df %>% 
  dplyr::mutate(.,lemma=rownames(.)) %>% 
  tidyr::gather(., key=lemma,value=value,na.rm=TRUE) 
food_LS_mat_long <- food_LS_mat_long[food_LS_mat_long$value!=0,]
colnames(food_LS_mat_long)[2] <- c("def")
write.table(food_LS_mat_long, "export/food_LS_defs.csv",sep = ",")

#Nodes and edges for visual inspection
def_LS_dist <- dist(food_LS_mat.df,"binary")
def_dist <- as.matrix(def_LS_dist) 
def_dist[which(def_dist == 0)] <- NA
links <- melt(def_dist, na.rm = TRUE)
rownames(links) <- NULL 
colnames(links) <- c("from","to","weight")
nodes <- unique(links[1])
net <- graph.data.frame(links,directed=T)
# The graph is irreadable (each node has a link to another one)
nodes_to_excl <- E(net)[E(net)$weight < mean(E(net)$weight)]
net_red <- delete_edges(net,nodes_to_excl)
l_circle <- layout.circle(net_red)
l_kamada <- layout.kamada.kawai(net_red)
l_spring <- layout.spring(net_red)
l_nicely <- layout_nicely(net_red)
l_frucht <- layout.fruchterman.reingold(net_red,repulserad=vcount(net_red)^1.5)
plot(net_red,
     edge.arrow.size=.4,vertex.label.cex=1.5,vertex.shape="none",
     edge.arrow.mode=0, edge.width = 0.1
     ,layout=l_frucht
     #,layout=l_nicely
     #,mark.groups=c(2,4)
     ,main = "Similarity of food-words accoding to the LS definitions"
     )

def_LS_dist <- dist(food_LS_mat.df)
def_LS_bin_dist <- dist(food_LS_mat.df.bin)

def_LS_cl <- hclust(def_LS_dist,
                    method ="ward.D")
def_LS_bin_cl <- hclust(def_LS_bin_dist,
                    method ="ward.D")

par(mfrow=c(1,1))
plot(def_LS_cl,main = "all data")
plot(def_LS_bin_cl,
     main="Word Clusters according to LS definitions",
     sub="(binary values)",
     xlab = "words", ylab="clusters")


#Finding nearest neighbours for a Latin term
syns <- nearest.neighbours(lemma_VObj_proj,
                           term = "alimentum", n=15, dist.matrix = T)
plot(syns)

#Find fuzzy
## Find which are not in the matrix
in_m <- intersect(food_LS$lemma,rownames(lemma_VObj_proj))
not_in_m <- food_LS$lemma[which( !food_LS$lemma %in% 
                                   in_m ) ]
## Generate comparison matrix
fuzzy_matrix <- adist(not_in_m,rownames(lemma_VObj_proj))
replace_cands <- sapply(seq(nrow(fuzzy_matrix)) , function(rowN) {
  rownames(lemma_VObj_proj)[which(fuzzy_matrix[rowN,] == min(fuzzy_matrix[rowN,]))]
})
names(replace_cands) <- not_in_m
replace_cands
accepted_cands <- c("cibaria1","cibaria2", "epule", "pastus", "uesco", "uisceratio")
no_cands <- c("alimentarius", "coliphia", "edendum", "pransus")
# No occurrences of no_cands in the texts:
## alimentarius: tylko alimentaria (1), alimentarie (1) +++
## coliphia: colyphium (1) - brak w MATRYCY! +++ 
## edendum: no occ (dubious in LS) ---
## pransus: (78) lemmatised as prandeo (807) +++

#accepted_cands <- c(accepted_cands, c("alimentaria", "colyphium","prandeo"))
accepted_cands <- c(in_m,accepted_cands)

# Wygenerować listę sąsiadów dla każdego wyrazu z listy LS: obrobić listę, by zawierała lemmaty z matrycy
food_LS_nearest <- lapply(accepted_cands, function(lemma) {
  nearest.neighbours(lemma_VObj_proj,
                         term = lemma, n=15, dist.matrix = T)
})
names(food_LS_nearest) <- accepted_cands
lapply(names(food_LS_nearest),function(lemma) {
  plot(food_LS_nearest[[lemma]],
       main=paste("Word similarity: ",lemma,""),
       sub="Patrologia Latina")
})

#Get list of the neigbouring lemma
food_PL_nearest <- data.frame(word = as.vector(
                                (sapply(names(food_LS_nearest),function(lemma) {
  neighbouring <- setdiff ( unlist( dimnames( food_LS_nearest[[lemma]] ) [1] ), lemma )
}))
)
)
food_PL_nearest_tbl <- dplyr::count(food_PL_nearest,word,sort = TRUE)

#Helper function: which word has some defined neighbourgh?
lapply(names(food_LS_nearest),function(x) which(unlist(dimnames(food_LS_nearest[[x]])[1]) == "elinguo"))
#Rysuje wykresy similarity dla wyrazów z PL
food_PL_nearest_sim <- lapply(as.vector(food_PL_nearest_tbl$word[1:30]), function(lemma) {
  nearest.neighbours(lemma_VObj_proj,
                     term = lemma, n=15, dist.matrix = T)
})

names(food_PL_nearest_sim) <- food_PL_nearest_tbl$word[1:30]

lapply( names(food_PL_nearest_sim), function(lemma) {
  plot(food_PL_nearest_sim[[lemma]],
       main=paste("Word similarity of 30 most frequent similar: ",lemma,""),
       sub="Patrologia Latina")
})

dplyr::top_n(food_PL_nearest_tbl,30)


#Measure intersection between LS and PL lists
intersect(accepted_cands,unlist(food_PL_nearest_tbl[,1])) 
write.table(food_PL_nearest_tbl,file = "export/food_frequent_neigbours.csv",
            col.names = TRUE, row.names = FALSE)
food_PL_neighb_ann <- read.csv(file="export/food_frequent_neigbours_ann.csv", 
                               header = T, sep = ",")
food_PL_neighb_ann <- food_PL_neighb_ann[,c("word", "n", "sense", "cat","subcat", "pos")]
#List of classified terms
count(food_PL_neighb_ann[,1:30],cat, subcat)
select(food_PL_neighb_ann[1:30,], cat, subcat)  %>% count(.,cat, sort = T)
## Food subcategories
select(food_PL_neighb_ann[1:30,], cat, subcat)  %>%
  filter(.,cat == "food") %>% 
  count(.,subcat, sort = T)
## Search for food of some category and subcategory
food_PL_neighb_ann[1:30,]  %>%  filter(cat=="food")  %>% 
  filter(subcat=="typeof")  %>% select(.,word)

## Clustering 30 most frequent words
freq_words <- food_PL_nearest_tbl$word[1:30]
freq_words_sim <- matrix(nrow = length(freq_words), ncol = length(freq_words))

freq_words_sim <- as.data.frame(freq_words_sim)

freq_words_sim <- lapply(freq_words, function(lemma1) {
  lapply(freq_words, function(lemma2) {
    pair.distances(w1 = lemma1, w2 = lemma2, M = lemma_VObj_proj)
  })
})

freq_words_sim <- as.data.frame(do.call(rbind, freq_words_sim), row.names = as.vector(freq_words))
colnames(freq_words_sim) <- as.vector(freq_words)
rownames(freq_words_sim) <- as.vector(freq_words)

freq_words.dist <- dist(freq_words_sim)
freq_words.PCA <- prcomp(freq_words.dist)
freq_word.clust <- hclust(freq_words.dist)
### PCA plot
plot(freq_words.PCA$x[,1],freq_words.PCA$x[,2], pch=NA,
     main="Clusters of Food Words", sub = "Patrologia Latin")
text(freq_words.PCA$x[,1],freq_words.PCA$x[,2], 
     labels=rownames(freq_words.PCA$x), pch=5, cex=1.2)
plot(freq_word.clust,
     main="Hierarchical Clusters of Food Words", sub = "Patrologia Latin")
library(dendextend)
freq_word.clust_dend <- as.dendrogram(freq_word.clust)
freq_word.clust_dend.cols <- sapply(labels(freq_word.clust_dend), function(label) { # Find category of each word
  food_PL_neighb_ann[food_PL_neighb_ann$word == label,"cat"]
} )
freq_word.clust_dend.cols <- freq_word.clust_dend.cols[drop=TRUE]
freq_word.clust_dend %>% 
  set("labels_col",as.numeric(freq_word.clust_dend.cols)) %>% 
  set("labels_cex", 1.7) %>% 
  hang.dendrogram %>% 
  plot(main="Sense Clusters of 30 Most Frequent Food Words",
       sub="Patrologia Latina"
       , horiz = FALSE) 
  legend("topright",
        legend = unique(freq_word.clust_dend.cols), 
         fill=as.numeric(unique(freq_word.clust_dend.cols)), cex = 1.3, y.intersp = 0.8)

#Frequency of neighbours
top30neighb_food <- as.data.frame(x=freq_words) %>%
  `colnames<-` ("lemma") %>%
  dplyr::mutate(.,freq=sapply(.[,1], function(word) {
  rcqp::cqi_id2freq("PATROLOGIA.lemma",
                    rcqp::cqi_str2id("PATROLOGIA.lemma",as.character(word) ))
})) %>% arrange(.,desc(freq)) 


#Export for gephi
tmp <- as.matrix(freq_words.dist)
tmp[tmp == 0] <- NA
tmp <- reshape2::melt(as.matrix(tmp),na.rm=TRUE)
colnames(tmp) <- c("Source", "Target", "Weight")
write.csv(tmp, file = "export/food_frequ_sim.csv")
rm(tmp)
#Next
#Określić relacje między znalezionymi wyrazami (czy da się sprawdzić ich pokrewieństwo: k-clusters, hclust?)
# Ludzkie i zwierzęce
# Czasowniki i rzeczowniki

library(MASS)
mds <- isoMDS(nn_l_p.dist, p=2)
pos <- sapply(dimnames(mds$points)[[1]],function(x) gsub(pattern = ".+_","",x))
pos_fac <- as.factor(pos) #Pos jako faktor służący nadaniu kolorów
plot(mds$points,pch=19, col=pos_fac)
text(mds$points, labels=names(nn_l_p.terms), pos=3,cex = 1.8)
legend(x="topright", legend = levels(pos_fac), col=pos_fac, pch=19)