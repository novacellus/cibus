library(data.table)
library(dplyr)
patro <- fread(input = "../wordspace/tables/patrologia.full.meta.tbl")
colnames(patro) <- c("freq","w1","id1","vol1","w2","id2","vol2")
patrologia.meta <- c("2", "4", "1", "19",
                     "4", "5", "20", "63",
                     "5", "6", "64", "73",
                     "6", "7", "74", "88",
                     "7", "8", "89", "96",
                     "8", "9", "97", "130",
                     "9", "10", "131", "136",
                     "10", "11", "137", "150",
                     "11", "12", "151", "174",
                     "12", "12", "175", "205",
                     "12", "13", "206", "217") #PL dating (teste Wikipedia)
names(patrologia.meta) <- rep(c("c.from","c.to","v.from","v.to"),11)

patrologia.meta <- data.frame(c.from=patrologia.meta[seq(1,length(patrologia.meta),4)],
                               c.to=patrologia.meta[seq(2,length(patrologia.meta),4)],
                               v.from=patrologia.meta[seq(3,length(patrologia.meta),4)],
                               v.to=patrologia.meta[seq(4,length(patrologia.meta),4)])
patrologia.time <- data.table(
  vol1=as.numeric(unlist(apply(patrologia.meta,1,function(x) seq(from = as.numeric(x[["v.from"]]),as.numeric(x[["v.to"]]), 1 )))),
  from=as.numeric(unlist(apply(patrologia.meta,1,function(x) rep(x[["c.from"]], times=(as.numeric(x[["v.to"]])-as.numeric(x[["v.from"]])+1) )))),
  to=as.numeric(unlist(apply(patrologia.meta,1,function(x) rep(x[["c.to"]],times=(as.numeric(x[["v.to"]])-as.numeric(x[["v.from"]])+1))) ))
  )

patro.time <- left_join(patro,patrologia.time,by="vol1") #Ab. 17s
View(patro.time[1:10,])
#Create subtables for similarity tables
patro.time[,sum(freq),by=c("to")] # Display word counts
# By century
patro.time.4 <- patro.time[to <= 4, sum(freq),by=c("w1","w2")]
patro.time.5 <- patro.time[to == 5, sum(freq),by=c("w1","w2")]
patro.time.6 <- patro.time[to == 6, sum(freq),by=c("w1","w2")]
patro.time.7 <- patro.time[to == 7, sum(freq),by=c("w1","w2")]
patro.time.8 <- patro.time[to == 8, sum(freq),by=c("w1","w2")]
patro.time.9 <- patro.time[to == 9, sum(freq),by=c("w1","w2")]
patro.time.10 <- patro.time[to == 10, sum(freq),by=c("w1","w2")]
patro.time.11 <- patro.time[to == 11, sum(freq),by=c("w1","w2")]
patro.time.12 <- patro.time[to == 12, sum(freq),by=c("w1","w2")]
patro.time.13 <- patro.time[to == 13, sum(freq),by=c("w1","w2")]

patro.time.tbls <- list(patro.time.4=patro.time.4,patro.time.5=patro.time.5,patro.time.6=patro.time.6,
                       patro.time.7=patro.time.7,patro.time.8=patro.time.8,patro.time.9=patro.time.9,
                       patro.time.10=patro.time.10,patro.time.11=patro.time.11,patro.time.12=patro.time.12,
                       patro.time.13=patro.time.13)
#Add names
lapply(patro.time.tbls,function(tbl) {
  colnames(tbl) <- c("w1","w2","freq")}
)

names(patro.time.tbls) <- c("patro.time.4","patro.time.5","patro.time.6","patro.time.7",
                           "patro.time.8","patro.time.9","patro.time.10","patro.time.11","patro.time.12","patro.time.13")
# Create similarity matrices for collocation tables
library(MASS)
library(wordspace)
#Not working
patro.time.sims_VObj <- lapply(patro.time.tbls,
                               function(lista) {
                                 dsm(target=lista$w1, feature=lista$w2,
                                     score=lista$freq, raw.freq=TRUE, sort=TRUE,verbose = T)
                               }
                               )
lapply(patro.time.sims_VObj, dim)

#3. Reducing matrix
patro.time.sims_VObj <- lapply(patro.time.sims_VObj, function(tbl) {
  subset(tbl, nnzero >= 3, nnzero >= 3, recursive=TRUE)
})

#b. Weighthing of Matrix
#Parameters:
patro.time.sims_VObj_weight <- lapply(patro.time.sims_VObj, function(matr) {
  dsm.score(matr, score="simple-ll", transform="log", normalize=TRUE)
            })

#c. Reduce by subspace projection (ca. 30')
patro_VObj_proj  <- lapply(patro.time.sims_VObj_weight, function(matr) {
  dsm.projection(matr, method="rsvd", n=300, oversampling=4,)
})
  
  
save(patro_VObj_proj,file = "patro_VObj_proj")

#load("patro_VObj_proj")

# Analysis
cibus.time.neig <-lapply(patro_VObj_proj,function(mat) nearest.neighbours(mat,"cibus")) %>% 
  lapply(names) %>% unlist #Similar lemmas
