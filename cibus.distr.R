library(data.table)
library(devtools)
library(rcqp)

load_all("../liber/lex/lexR/")

breaks.n <- 100

corp.size <- rcqp::cqi_attribute_size(paste("PATROLOGIA",".",
                                            "lemma",sep=""))
thresh <- seq (from = 1, to=corp.size,
               by = corp.size / breaks.n ) %>%
  c(., corp.size) #Calculates breaks at which data will be cut

breaks.count <- round(corp.size / breaks.n) #How much tokens each break contains

cibus.pos <- find_word_pos(corpus = "PATROLOGIA",attr = "lemma",what.id=find_word_id(corpus = "PATROLOGIA",attr = "lemma",what = "cibus"))
cibus.pos.binned <- cut(x = cibus.pos, breaks = thresh, labels = 1:breaks.n)

cibus.binned <- data.table(pos=cibus.pos,bin=cibus.pos.binned)
cibus.by.bin <- cibus.binned[,nrow(.SD),by=bin] 
colnames(cibus.by.bin) <- c("bin","count")

cibus.pos.2struc <- cqi_cpos2struc(attribute = "PATROLOGIA.text_vol",cpos = cibus.pos)

cibus.binned.meta <- cibus.binned %>% 
  mutate(vol=cqi_struc2str(attribute = "PATROLOGIA.text_vol",cibus.pos.2struc),
         aut=cqi_struc2str(attribute = "PATROLOGIA.text_auteur",cibus.pos.2struc),
         tit=cqi_struc2str(attribute = "PATROLOGIA.text_titre",cibus.pos.2struc)
  )
setkey(cibus.binned.meta,bin,vol)

cibus.by.bin <- dplyr::mutate(cibus.by.bin,perM=(count/breaks.count)*1000000) #Words per million
cibus.labels <- sapply (cibus.by.bin$bin, function(binn) {
    paste(unique(cibus.binned.meta[as.character(binn),vol]),collapse = "_")#Get the list of vols/aut etc. in a respective bin
})

cibus.by.max <- cibus.by.bin[order(cibus.by.bin$perM)][91:100] #Bins of 10 max
labels.max <- cibus.labels[cibus.by.max$bin]
aut.max <- unique(cibus.binned.meta[as.vector(cibus.by.max[10:1][2]$bin),aut])
tit.min <- unique(cibus.binned.meta[as.vector(cibus.by.max[10:1][2]$bin),nrow(.SD),by=tit])
aut.tit.max <- unique(cibus.binned.meta[as.vector(cibus.by.max[10:1][2]$bin),nrow(.SD),by=list(aut,tit)]) %>% 
  mutate(aut_tit=paste(aut,tit,sep="_"))
aut.tit.max[order(aut.tit.max$V1,decreasing = T)][1:5]

### Max: 73-74-75
#1:                           Theodoretus Cyrensis      Philotheus, sive Theophiles 110
#2:       Heraclides Alexandrinus; Auctor incertus             Heraclidis Paradisus  53
#3:                      Palladius Helenopolitanus                Historia Lausiaca  28
#4:             Auctor incertus; Pelagius diaconus                   Verba seniorum  27
#5:           Joannes Damascenus; Billius, Jacobus     Vita SS. Barlaam et Josaphat  19

### Max: 176-177-178
#1: Auctor incertus (Hugo de S. Victore?)               Miscellanea 87
#2:                     Petrus Abaelardus                 Epistolae 82
#3:    Auctor incertus (Hugo de Folieto?) De bestiis et aliis rebus 71
#4:                       Auctor incertus                  Sermones 54
#5:                     Petrus Abaelardus                    Ethica 22

cibus.by.min <- cibus.by.bin[order(cibus.by.bin$perM)][1:10] #Bins of 10 min
labels.min <- cibus.labels[cibus.by.min$bin]
aut.min <- unique(cibus.binned.meta[as.vector(cibus.by.min[1:10][2]$bin),aut])
tit.min <- unique(cibus.binned.meta[as.vector(cibus.by.min[1:10][2]$bin),nrow(.SD),by=tit])
aut.tit.min <- unique(cibus.binned.meta[as.vector(cibus.by.min[1:10][2]$bin),nrow(.SD),by=list(aut,tit)]) %>% 
  mutate(aut_tit=paste(aut,tit,sep="_"))
aut.tit.min[order(aut.tit.min$V1,decreasing = T)]
### Min: 214-215
#1: Auctor incertus   Gesta Innocentii III  1   Auctor incertus_Gesta Innocentii III
#2: Innocentius III Regesta sive epistolae 27 Innocentius III_Regesta sive epistolae

### Min: 126-130
#1:                          Isidorus Mercator                  Collectio decretalium 15
#2:                               Joannes VIII                   Epistolae et decreta  9
#3:                        Hincmarus Rhemensis Opuscula in causa Hincmari Laudunensis  7
#4: Auctor incertus (Liutprandus Cremonensis?)          De pontificum Romanorum vitis  5
#5: Anastasius bibliothecarius; Auctores varii                            Collectanea  4
#6:                            Auctor incertus                             De computo  4
#7:                   Wolfhardus Hasenrietanus                      Vita S. Walpurgis  3
#8:                        Hincmarus Rhemensis                              Epistolae  1
#9:                 Anastasius bibliothecarius    Interpretatio Synodi VIII generalis  1
#10:                 Anastasius bibliothecarius     Interpretatio Synodi VII generalis  1
#11:  Auctor incertus (Monachus Ratisbonensis?)                         Benedictio Dei  1
#12:                            Auctor incertus             Invectiva pro Formoso papa  1
#13:                        Hebernus Turonensis                    Miracula B. Martini  1
#14:                                Stephanus V     Epistolae, diplomata et privilegia  1

cibus.by.bin %>% ggplot() + aes(x=bin,y=perM) + 
  geom_point(size=5) +  guides(size = "legend", colour = "none") +
  labs(title = "CIBUS Distribution\nPatrologia Latina") +
  xlab("Corpus Position") + 
  theme(axis.text.x = element_text(angle = 270)) +
  ylab("Occurrences per million words") + 
  geom_point(data = cibus.by.max, 
             aes(x=bin,y=perM,label=labels.max, color="red", size=6)) +
  geom_point(data = cibus.by.min, 
             aes(x=bin,y=perM,label=labels.min, color="green", size=6)) +
  geom_text(data = cibus.by.max, 
            aes(x=bin,y=perM,label=labels.max, color="red")) + 
  geom_text(data = cibus.by.min, 
            aes(x=bin,y=perM,label=labels.min, color="green"))

patrologia.meta.full <- fread("import/patrologia.metadane.full",header = F,sep = "\t")
colnames(patrologia.meta.full) <- c("n","id","vol","aut", "tit", "per")
patrologia.meta.full %>% arrange(patrologia.meta.full,id)

### Per text
patro.cibus.texts <- fread("import/patro.cibus.PerText.txt")
colnames(patro.cibus.texts) <- c("id","n","hits","perMil")
patro.cibus.texts <- left_join(patro.cibus.texts,patrologia.meta.full[,!"n",with=F],by="id")
top <- patro.cibus.texts[n<5000 & hits > 5][order(-perMil)][1:20] #First 20 texts with the largest number of CIBUS occurrences
write.table(top[,list(hits,perMil,n,aut,tit)],file = "export/tmp.txt")

