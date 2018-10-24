# annotations

## I didn't think Arrystar annotation was great, but I included it. 
```{r eval=FALSE}
load("agilannotAll.rda")
```


I hope you use this, because I spent alot of time on it and it would make me happy if it helped someone else.
I used the code from this paper [Re-Annotator: Annotation Pipeline for Microarrays](https://www.biorxiv.org/content/early/2015/05/21/019596) to align the probe sequences to the exome. I did this for several databases and Hg19 and 38 and then combined them. 

## The end list is the first hit from this process
```{r, eval=FALSE}
load("firstHit.rda")
```

## I also included the files that have all the hits:
```{r, eval=FALSE}
# matches for each database/build
load("ReAnnoLists.rda")

# datafram with all hits in one field (separated by ";") after compiling from different databases/builds
load("MergedReAnno.rda")

# list with dataframe of separate hits from the merged table
load("HitList.rda")

# dataframe of the first hit for each probe
load("firstHit.rda")
```

code for compiling the output to help explain the files.
```{r eval=FALSE}
# import and add variables for build, database, ID types
gcComp28<-read.delim("GCcomp28_probesAgilAll_exome_readAnnotation.txt", row.names = 1)
gcComp28$ProbeID<-rownames(gcComp28)
gcComp28$build<-"Hg38"; gcComp28$Db<-"wgGenCodeCompV28"
gcComp28$TxID<-"ENST"
gcComp28$AltID<-"gene_id"

kg38<-read.delim("knowngene38_probesAgilAll_exome_readAnnotation.txt", row.names = 1)
kg38$ProbeID<-rownames(kg38)
kg38$build<-"Hg38"; kg38$Db<-"knowngenes"
kg38$TxID<-"UCSC"
kg38$AltID<-"ENST"

refseq38<-read.delim("RefSeq38_probesAgilAll_exome_readAnnotation.txt", row.names = 1)
refseq38$ProbeID<-rownames(refseq38)
refseq38$build<-"Hg38"; refseq38$Db<-"RefSeq"
refseq38$TxID<-"RefSeq"
refseq38$AltID<-"gene_id"

kg19<-read.delim("knowngene19_probesAgilAll_exome_readAnnotation.txt", row.names = 1)
kg19$ProbeID<-rownames(kg19)
kg19$build<-"Hg19"; kg19$Db<-"knowngenes"
kg19$TxID<-"UCSC"
kg19$AltID<-"UCSC"

refseq19<-read.delim("RefSeq19_probesAgilAll_exome_readAnnotation.txt", row.names = 1)
refseq19$ProbeID<-rownames(refseq19)
refseq19$build<-"Hg19"; refseq19$Db<-"RefSeq"
refseq19$TxID<-"RefSeq"
refseq19$AltID<-"gene_id"

ens19<-read.delim("ensGene19_probesAgilAll_exome_readAnnotation.txt", row.names = 1)
ens19$ProbeID<-rownames(ens19)
ens19$build<-"Hg19"; ens19$Db<-"ensGene"
ens19$TxID<-"ENST"
ens19$AltID<-"ENST"

ReAnnoList_Hg19<-list(ENS=ens19,
                      UCSC=kg19,
                      NCBI=refseq19)

ReAnnoList_Hg38<-list(GENCODE=gcComp28,
                      UCSC=kg38,
                      NCBI=refseq38)

save(ReAnnoList_Hg19, ReAnnoList_Hg38, file = "ReAnnoLists.rda")
lapply(ReAnnoList_Hg19, nrow)
# $ENS
# [1] 40568
# 
# $UCSC
# [1] 31251
# 
# $NCBI
# [1] 30819

lapply(ReAnnoList_Hg38, nrow)
# $GENCODE
# [1] 39148
# 
# $UCSC
# [1] 39196
# 
# $NCBI
# [1] 34470
# Build composite data 

reanno<-ReAnnoList_Hg38$GENCODE 
reanno<-rbind(reanno,ReAnnoList_Hg19$ENS[!rownames(ReAnnoList_Hg19$ENS)%in%rownames(reanno),])
nrow(reanno) #41159

reanno<-rbind(reanno,ReAnnoList_Hg38$UCSC[!rownames(ReAnnoList_Hg38$UCSC)%in%rownames(reanno),])
nrow(reanno) #41169

reanno<-rbind(reanno,ReAnnoList_Hg38$NCBI[!rownames(ReAnnoList_Hg38$NCBI)%in%rownames(reanno),])
nrow(reanno) #43352

reanno<-rbind(reanno,ReAnnoList_Hg19$UCSC[!rownames(ReAnnoList_Hg19$UCSC)%in%rownames(reanno),])
nrow(reanno) #43959

reanno<-rbind(reanno,ReAnnoList_Hg19$NCBI[!rownames(ReAnnoList_Hg19$NCBI)%in%rownames(reanno),])
nrow(reanno) #43962

save(reanno, file = "MergedReAnno.rda")

# reanno has one HIT field with multiple hits separated by ";" 
# HitList is List with dataframe of hits for each probe

HitList<-strsplit(reanno$HIT,"\\;")
names(HitList)<-reanno$ProbeID

HitList<-lapply(HitList,function(x){
        x<-gsub("\\,","\\|",x)
        x<-as.data.frame(do.call(rbind,strsplit(x,"\\|")),stringsAsFactors=FALSE)
        names(x)<-c("chr","st","ranges","TxID","AltID","region")
        return(x)
})

#removes numbers in front of IDs
HitList<-lapply(HitList, function(x){
        x[,"TxID"]<-gsub("[0-9]+_","",x[,"TxID"])
        return(x)
})

save(HitList,file = "HitList.rda")

#Make one table with just first hit
firstHit<-lapply(HitList, function(x){x[1,]})
firstHit<-do.call(rbind,firstHit)
save(firstHit, file = "firstHit.rda")
```
