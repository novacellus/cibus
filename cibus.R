#Data
cibus.colloc <- read.csv("import/collocs/cibus.txt",header = TRUE, sep = "\t")
colnames(cibus.colloc) <- c("rank", "lemma", "freq", "exp", "obs", "no_texts", "dice")
cibus.colloc <- cibus.colloc %>%
  mutate(., pos= sapply(lemma, function(word) {
    cqi_cpos2str("PATROLOGIA.pos",
      cqi_id2cpos("PATROLOGIA.lemma",
        cqi_str2id("PATROLOGIA.lemma",as.vector(word) ))[1])
}))


cibus.neighs1000 <- nearest.neighbours(M = lemma_VObj_proj,"cibus",n = 1000)
cibus.neighs1000.df <- as.data.frame(cibus.neighs1000) %>% 
  mutate(.,lemma=rownames(cibus.neighs1000.df))
#Overlap between collocates and neighbours
overlap <- sapply(1:1000,function(x) {
  length(intersect(cibus.colloc$lemma[1:x],names(cibus.neighs1000[1:x])))
})
overlap.df <- data.frame(count=1:1000, overlap=overlap)
ggplot(overlap.df,aes(x=count,y=overlap,colour="real overlap")) + 
  geom_point() + 
  geom_line(aes(x=1:1000,y=1:1000,colour="theoretical overlap")) + 
  ggtitle("Similar Words vs. Collocates\ncibus in Patrologia Latina") + 
  theme(legend.position="bottom")

not_in_neighs <- setdiff(cibus.colloc$lemma[1:1000], names(cibus.neighs1000))
not_in_collocs <- setdiff(names(cibus.neighs1000),cibus.colloc$lemma[1:1000])
in_both <- intersect(names(cibus.neighs1000),cibus.colloc$lemma[1:1000])

#Creates links for net visuals: rank column is word's rank with 100 attributed to the best
sim_net <- data.frame(source=
                        c(as.character(cibus.colloc$lemma[1:50]),
                          names(cibus.neighs1000)[1:50]),
                      target=c(rep("collocation",100), rep("neighbourgh",100)),
                      rank=c(100:1,100:1) )
overlap.net <- graph.data.frame(sim_net,directed = F)

l_circle <- layout.circle(overlap.net)

#cols <- c("red","blue")
#E(overlap.net)$color <- cols[E(overlap.net)$target]
plot(overlap.net, 
     vertex.label.cex=1.0,vertex.shape="none"
     ,layout=layout.circle
     )

#Another attempt
nodes.id <- c(
  sapply (as.character(cibus.colloc$lemma[1:100]), function(x) paste(x,"_colloc",sep="")),
  sapply(names(cibus.neighs1000)[1:100], function(x) paste(x,"_neighb",sep=""))
         )
nodes.lemma <- c(as.character(cibus.colloc$lemma[1:100]), names(cibus.neighs1000)[1:100])
nodes.rank <- c(100:1,100:1)
nodes.type <- c(rep("collocation",100), rep("neighbourgh",100))
nodes <- data.frame(id=nodes.id,lemma=nodes.lemma,
                    rank=nodes.rank, type=nodes.type)

overl.nodes <- intersect(nodes$lemma[1:100], nodes$lemma[101:200])
links <- data.frame(
  from = as.vector(sapply(overl.nodes, function(word1) {
  nodes[1:100,][nodes[1:100,]$lemma == word1, "id"]
})),
  to = as.vector(sapply(overl.nodes, function(word1) {
    nodes[101:200,][nodes[101:200,]$lemma == word1, "id"]
  })),row.names = NULL
)

overlap.net2 <- graph.data.frame(d = links, vertices = nodes,directed = T)

#Write to file
write.csv(links,file = "export/cibus_overlap.links.csv")
write.csv(nodes,file = "export/cibus_overlap.nodes.csv")

#Excluding nodes
#Only rank > 50: del.nodes <- V(overlap.net2)[V(overlap.net2)$rank<51]
#Only having links
del.nodes <- V(overlap.net2)[!V(overlap.net2)$lemma %in% intersect(nodes$lemma[1:100], nodes$lemma[101:200])]
overlap.net2.red <- delete_vertices(overlap.net2,del.nodes)
#Parallel layout
l <- matrix(c(c(rep(1,32), rep(10,32)), c(1:32,1:32) ), ncol=2 )
l1 <- layout.fruchterman.reingold(overlap.net2.red)
l2 <- layout.circle(overlap.net2.red)
cols=c("red","green")
V(overlap.net2.red)$color <- cols[as.factor(V(overlap.net2.red)$type)]
plot(overlap.net2.red,layout=l,vertex.shape="none"
     ,edge.arrow.mode=0,vertex.label.color=V(overlap.net2.red)$color
     ,vertex.label=c(sort(V(overlap.net2.red)$lemma[1:(vcount(overlap.net2.red)/2)],decreasing = T),
                     sort(V(overlap.net2.red)$lemma[(vcount(overlap.net2.red)/2+1):vcount(overlap.net2.red)],decreasing = T)
     ),vertex.label.cex=1,edge.width=0
     ,main="Overlap between similar and neighbourgh words\nCIBUS"
     ,sub="First 100 hundred words"
     )

#Overlap between nouns only
t <- read.csv("export/cibus.neighs_ann.csv")
edg <- cbind(as.vector(cibus.colloc[cibus.colloc$pos=="SUB",]$lemma[1:100]),
             as.vector(t[t$pos=="N","lemma"][1:100]))
ovlps <- intersect(as.vector(cibus.colloc[cibus.colloc$pos=="SUB",]$lemma[1:100]),
                   as.vector(t[t$pos=="N","lemma"][1:100]))

tt.nod = data.frame(
  c(sapply(edg[,1], function(x) paste(x,"_col",sep = "")),
    sapply(edg[,2], function(x) paste(x,"_sim",sep = "")))
)

tt.edg <- data.frame(from=sapply(sort(ovlps), function(x) paste(x,"_col",sep = "")),
                     to=sapply(sort(ovlps), function(x) paste(x,"_sim",sep = ""))
)

tt.net <- graph.data.frame(d = tt.edg, vertices = tt.nod)
l <- matrix(c(c(rep(1,100), rep(5,100)), c(100:1,100:1) ), ncol=2 )
cols=c("red","green")
V(tt.net)$color <- c(rep("red",100),rep("green",100))
plot(tt.net,layout=l,edge.arrow.mode=0,vertex.shape="none",vertex.label.color=V(tt.net)$color,
     main="Overlap between noun collocations and similar words")


cibus.top <- cibus.colloc[1:30,]
top15_neigh <- nearest.neighbours(lemma_VObj_proj,
                   term = "cibus", n=30, dist.matrix = F)

#Clustering on neighbours and collocations
## CLustering on neighbours
howmany <- 30
cibus.neighs.100 <- names(cibus.neighs1000[1:howmany])
cibus.sim.df <- matrix(nrow = length(cibus.neighs.100),
                       ncol = length(cibus.neighs.100)) %>% 
  as.data.frame
cibus.sim.df <- lapply(cibus.neighs.100, function(lemma1) {
  lapply(cibus.neighs.100, function(lemma2) {
    pair.distances(w1 = lemma1, w2 = lemma2, M = lemma_VObj_proj)
  })
})
cibus.sim.df <- as.data.frame(do.call(rbind, cibus.sim.df),
                              row.names = as.vector(cibus.neighs.100))
cibus.sim.df <- sapply(cibus.sim.df,as.numeric) #Makes sure all values are numerical
colnames(cibus.sim.df) <- as.vector(cibus.neighs.100)
rownames(cibus.sim.df) <- as.vector(cibus.neighs.100)

cibus.sim.dist <- dist(cibus.sim.df)
cibus.sim.PCA <- prcomp(x = cibus.sim.df,na.rm=FALSE)
cibus.sim.clust <- hclust(cibus.sim.dist)
### PCA plot
plot(cibus.sim.PCA$x[,1],cibus.sim.PCA$x[,2], pch=NA,
     main="Clusters of CIBUS related words", sub = "Patrologia Latina")
text(cibus.sim.PCA$x[,1],cibus.sim.PCA$x[,2], 
     labels=rownames(cibus.sim.PCA$x), pch=5, cex=1.2)
plot(cibus.sim.clust,
     main="Hierarchical Clusters of CIBUS Related Words", sub = "Patrologia Latina")
library(dendextend)
cibus.sim.clust_dend <- as.dendrogram(cibus.sim.clust)
cibus.sim.clust_dend %>% 
  hang.dendrogram %>% 
  plot(main="Sense Clusters of 30 CIBUS related Words",
       sub="Patrologia Latina"
       , horiz = TRUE, axes=FALSE) 
#Flattene-out
par(mfrow=c(1,1))
cibus.sim.clust_dend %>% 
  collapse_branch(.,tol=0.9,lower = TRUE) %>% 
  #ladderize %>% 
  hang.dendrogram %>% 
  set("labels_colors") %>% 
  plot(main="Sense Clusters of 30 CIBUS related Words",
       sub="Patrologia Latina"
       , horiz = TRUE, axes=FALSE) 

## CLustering on collocations
howmany <- 50
cibus.collocs.50 <- as.vector(cibus.colloc[1:howmany,"lemma"])
cibus.coll.sim.df <- matrix(nrow = length(cibus.collocs.50),
                       ncol = length(cibus.collocs.50)) %>% 
  as.data.frame
cibus.coll.sim.df <- lapply(cibus.collocs.50, function(lemma1) {
  lapply(cibus.collocs.50, function(lemma2) {
    pair.distances(w1 = lemma1, w2 = lemma2, M = lemma_VObj_proj)
  })
})
cibus.coll.sim.df <- as.data.frame(do.call(rbind, cibus.coll.sim.df),
                              row.names = as.vector(cibus.collocs.50))
cibus.coll.sim.df <- sapply(cibus.coll.sim.df,as.numeric) #Makes sure all values are numerical
colnames(cibus.coll.sim.df) <- as.vector(cibus.collocs.50)
rownames(cibus.coll.sim.df) <- as.vector(cibus.collocs.50)

cibus.coll.sim.dist <- dist(cibus.coll.sim.df)
cibus.coll.sim.PCA <- prcomp(x = cibus.coll.sim.df,na.rm=FALSE)
cibus.coll.sim.clust <- hclust(cibus.coll.sim.dist)
### PCA plot
plot(cibus.coll.sim.PCA$x[,1],cibus.coll.sim.PCA$x[,2], pch=NA,
     main="Clusters of CIBUS related words", sub = "Patrologia Latina")
text(cibus.coll.sim.PCA$x[,1],cibus.coll.sim.PCA$x[,2], 
     labels=rownames(cibus.coll.sim.PCA$x), pch=5, cex=1.2)
plot(cibus.coll.sim.clust,
     main="Hierarchical Clusters of CIBUS Related Words", sub = "Patrologia Latina")
library(dendextend)
cibus.coll.sim.clust_dend <- as.dendrogram(cibus.coll.sim.clust)
cibus.coll.sim.clust_dend %>% 
  hang.dendrogram %>% 
  plot(main="Sense Clusters of 50 CIBUS Collocations",
       sub="Patrologia Latina"
       , horiz = TRUE, axes=FALSE) 
#Flatten-out
par(mfrow=c(1,1))
cibus.coll.sim.clust_dend %>% 
  collapse_branch(.,tol=0.9,lower = TRUE) %>% 
  ladderize %>% #Ladderizing pokazuje jako pierwsze węzły o najmniejszej liczbie węzłów podrzędnychO 
  hang.dendrogram %>% 
  set("labels_colors") %>%
  #prune(cibus.colloc[1:howmany,][cibus.colloc[1:howmany,]$pos != "VBE","lemma"]) %>% 
  plot(main="Sense Clusters of 50 CIBUS Collocations",
       sub="Patrologia Latina"
       , horiz = FALSE, axes=FALSE)
cibus.coll.sim.clust_dend %>% ladderize %>% 
  rect.dendrogram(k=2, border = 8, lty = 3, lwd = 3, 
                  horiz=FALSE,text=c("NOUNS AND ADJECTIVES", "VERBS") ) %>% 
  
  #VERBS ONLY
cibus.coll.sim.clust_dend %>% 
  collapse_branch(.,tol=0.9,lower = TRUE) %>% 
  ladderize %>% 
  hang.dendrogram %>% 
  set("labels_colors") %>%
  prune(cibus.colloc[1:howmany,][cibus.colloc[1:howmany,]$pos != "VBE","lemma"]) %>% 
  plot(main="Sense Clusters of 50 CIBUS Collocations\nVERBS",
       sub="Patrologia Latina"
       , horiz = TRUE, axes=FALSE, cex=2)
##NOUNS
cibus.coll.sim.clust_dend %>% 
  collapse_branch(.,tol=0.9,lower = TRUE) %>% 
  ladderize %>% 
  hang.dendrogram %>% 
  set("labels_colors") %>%
  prune(cibus.colloc[1:howmany,][cibus.colloc[1:howmany,]$pos == "VBE" | cibus.colloc[1:howmany,]$pos == "QLF","lemma"]) %>% 
  plot(main="Sense Clusters of 50 CIBUS Collocations\nNOUNS",
       sub="Patrologia Latina"
       , horiz = TRUE, axes=FALSE, cex=2)

SUB <- as.vector(cibus.colloc[1:howmany,][cibus.colloc[1:howmany,]$pos == "SUB","lemma"])
VBE <- as.vector(cibus.colloc[1:howmany,][cibus.colloc[1:howmany,]$pos == "VBE","lemma"])
ADJ <- as.vector(cibus.colloc[1:howmany,][cibus.colloc[1:howmany,]$pos == "QLF","lemma"])
words_to_exclude = c("solidum2","huiusmodi")

##ADJ
cibus.coll.sim.clust_dend %>% 
  collapse_branch(.,tol=0.9,lower = TRUE) %>% 
  ladderize %>% 
  hang.dendrogram %>% 
  set("labels_colors") %>%
  prune(c(VBE,SUB) ) %>% #WOrds to exclude from graph
  plot(main="Sense Clusters of 50 CIBUS Collocations\nADJECTIVES",
       sub="Patrologia Latina"
       , horiz = TRUE, axes=FALSE, cex=2)

#Comparison of collocations
panis.colloc <- read.csv("import/collocs/panis.txt",header = TRUE, sep = "\t")
colnames(panis.colloc) <- c("rank", "lemma", "freq", "exp", "obs", "no_texts", "dice")
panis.colloc <- panis.colloc %>%
  mutate(., pos= sapply(lemma, function(word) {
    cqi_cpos2str("PATROLOGIA.pos",
                 cqi_id2cpos("PATROLOGIA.lemma",
                             cqi_str2id("PATROLOGIA.lemma",as.vector(word) ))[1])
  }))
panis.lemmas <- as.vector(panis.colloc$dice)
names(panis.lemmas) <- panis.colloc$lemma
cibus.lemmas <- as.vector(cibus.colloc$dice)
names(cibus.lemmas) <- cibus.colloc$lemma

collocs.all <- sort(union(panis.colloc$lemma,cibus.colloc$lemma))
collocs.1 <- panis.lemmas[collocs.all] #Results for 1st colloc
collocs.2 <- cibus.lemmas[collocs.all] #Results for 2nd colloc
collocs.diff <- data.table(lemma=collocs.all,col_1=collocs.1, col_2=collocs.2)

collocs.diff[which(is.na(col_1),col_1),col_1 := -min(collocs.diff[col_1 > 0, col_1])] # Replace NAs with -(min effect other than 0) (waiting for better solution)
collocs.diff[which(is.na(col_2),col_2),col_2 := -min(collocs.diff[col_2 > 0, col_2])]

collocs.diff[, diff1 := col_1-col_2] #Differences for 'panis'
collocs.diff[, diff2 := col_2-col_1] #Differences for 'cibus'

collocs.1.diff <- collocs.diff[order(-diff1)][1:20] #Top 20 differences for 'panis'
collocs.2.diff <- collocs.diff[order(-diff2)][1:20] #Top 20 differences for 'cibus'
collocs.sim <- collocs.diff[col_1>0.005 & col_1 == col_2]

write.table(rbind(collocs.1.diff,collocs.2.diff),file = "export/panis.cibus.collocs.diff.csv")


cibus.id <- cqi_str2id("PATROLOGIA.lemma","cibus")
cibus.cpos <- cqi_id2cpos("PATROLOGIA.lemma",cibus.id)
cibus.struc <- cqi_cpos2struc("PATROLOGIA.text_id",cibus.cpos)
cibus.texts <- cqi_struc2str("PATROLOGIA.text_id",cibus.struc)
cibus.all <- data.table(id=cibus.texts) %>% 
  .[, .N, by=id] %>% left_join(.,patrologia.meta.full,by="id") %>% 
  .[,cibus.pmw := (N/n)*1000000]

panis.id <- cqi_str2id("PATROLOGIA.lemma","panis")
panis.cpos <- cqi_id2cpos("PATROLOGIA.lemma",panis.id)
panis.struc <- cqi_cpos2struc("PATROLOGIA.text_id",panis.cpos)
panis.texts <- cqi_struc2str("PATROLOGIA.text_id",panis.struc)
panis.all <- data.table(id=panis.texts) %>% 
  .[, .N, by=id] %>% left_join(.,patrologia.meta.full,by="id") %>% 
  .[,panis.pmw := (N/n)*1000000]

#Excluding some of cols
#panis.cibus.all <- full_join(panis.all[,.(n,panis.pmw,id)],cibus.all[,.(n,cibus.pmw,id)], by="id") %>% #TODO: exclude cols
#  full_join(patrologia.meta.full,.,by="id") %>% as.data.table
#Without excluding cols (better for further filtering)
panis.cibus.all <- full_join(panis.all,cibus.all, by="id") %>% #TODO: exclude cols
  full_join(patrologia.meta.full,.,by="id") %>% as.data.table

panis.cibus.all[which(is.na(panis.cibus.all[,cibus.pmw])),cibus.pmw := 0]
panis.cibus.all[which(is.na(panis.cibus.all[,panis.pmw])),panis.pmw := 0]
panis.cibus.all[,diff := (panis.pmw-cibus.pmw)] #Add diff (panis-cib)

cibus.diff <- panis.cibus.all[n>5000 & (n.x > 5 | n.y > 5)][order(diff)][1:20,.(aut,tit,panis.pmw,cibus.pmw)] #Small panis, big cibus (N.x - number of panis in the text)
panis.diff <- panis.cibus.all[n>5000 & (n.x > 5 | n.y > 5)][order(-diff)][1:20,.(aut,tit,panis.pmw,cibus.pmw)] #Big panis, small cibus
write.table(rbind(panis.diff,cibus.diff),"export/cibus.panis.txts.pmw.csv")

with(data = panis.cibus.all,cor(x=panis.pmw,y=cibus.pmw)) #Correlation of both frequencies
with(data = panis.cibus.all[(panis.pmw==0 & cibus.pmw!=0) | (panis.pmw!=0 & cibus.pmw==0) | (panis.pmw!=0 & cibus.pmw!=0)],
     cor(x=panis.pmw,y=cibus.pmw)) #Correlation of both frequencies
ggplot(panis.cibus.all,aes(x=panis.pmw,y=cibus.pmw)) + 
  geom_point()
par(mfrow=c(1,2))
#Including cases where both are zeros
ggplot(panis.cibus.all,aes(x=panis.pmw,y=cibus.pmw)) + 
  geom_point() + xlab("panis") + ylab("cibus") + 
  labs(title="Correlation of distribution\n0s included")
#Excluding cases where both are zeros
to_excl <- "(panis.pmw==0 & cibus.pmw!=) | (panis.pmw!=0 & cibus==0) | (panis.pmw!=0 & cibus!=0)"
ggplot(panis.cibus.all[(panis.pmw==0 & cibus.pmw!=0) | (panis.pmw!=0 & cibus.pmw==0) | (panis.pmw!=0 & cibus.pmw!=0)],
       aes(x=panis.pmw,y=cibus.pmw)) + 
  geom_point() + xlab("panis") + ylab("cibus") + 
  labs(title="Correlation of distribution\n0s excluded")


