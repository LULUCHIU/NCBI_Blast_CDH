library(biomaRt)
library(UniprotR)
library(vctrs)
library(stringr)
library(openxlsx)

# BiocManager::install("vctrs")
# BiocManager::install("UniprotR",force=TRUE)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

alist<-listAttributes(ensembl) # check all attributes

#> 1. Extract Extracellular Domain length sequence

familylist<-read.xlsx("/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/caherin_superfamily_homo_sapiens_bar.xlsx")

ENTREZID<-getBM(attributes=c('hgnc_symbol','entrezgene_id','uniprotswissprot'), 
                filters = 'entrezgene_id', 
                values = familylist$GeneID, 
                mart = ensembl)
ENTREZID<-ENTREZID[-which(ENTREZID$uniprotswissprot==""),]

addfix<-data.frame(hgnc_symbol=c("EPHB6","TGFB1"),entrezgene_id=c(2051,7040),uniprotswissprot=c("O15197",
                                                                                                "P01137"))
ENTREZID<-data.frame(rbind(ENTREZID,addfix))

CDHfamily<-data.frame()
for(jj in 1:length(ENTREZID$hgnc_symbol)){

id<-ENTREZID$uniprotswissprot[jj]
symbol<-ENTREZID$hgnc_symbol[jj]

domaininfo<-GetProteinAnnontate(id,"ft_domain")
ori_seq<-GetProteinAnnontate(id,"sequence")

if(domaininfo=="NA"){
  domainposi<-data.frame(symbol=symbol,Domain=c("NA"),start=c("NA"),end=c("NA"),
                         Name=c("NA"),Info=c("NA"),seq=c("NA"))
}else{

domain_list<-unlist(str_split(domaininfo,"; D"))

tmplist<-str_split(domain_list,";")

tmplist<-data.frame(do.call("rbind", tmplist))
if(dim(tmplist)[2]==2){
  tmplist<-data.frame(cbind(tmplist,data.frame(X3="")))
}
colnames(tmplist)<-c("Domain","Name","info")
domainposi<-data.frame(do.call("rbind",str_split(gsub("[A-Z]+ ","",tmplist$Domain),"\\..")))
domainposi<-cbind(domainposi,gsub(" \\/note=","",tmplist$Name))
domainposi<-cbind(domainposi,tmplist$info)
colnames(domainposi)<-c("start","end","Name","Info")
domainposi$seq<-NA
for(ii in 1:(dim(domainposi)[1])){
domainposi$seq[ii]<-substr(ori_seq, start = as.numeric(domainposi$start[ii]), 
                           stop = as.numeric(domainposi$end[ii]))
}
tmpc<-data.frame(symbol=rep(symbol,dim(domainposi)[1]),Domain=paste0("Domain ",c(1:dim(domainposi)[1])))
domainposi<-data.frame(cbind(tmpc,domainposi))
}
CDHfamily<-data.frame(rbind(CDHfamily,domainposi))
print(jj)
}

All<-merge(ENTREZID,CDHfamily,by.x="hgnc_symbol",by.y="symbol",all.y=T)
write.xlsx(All,"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/CDHfamily_domain_seq.xlsx")


All_clean<-All[-which(All$Domain=="NA"),]
All_clean$Domain<-gsub(" ","",All_clean$Domain)
All_clean$ID<-paste0(All_clean$hgnc_symbol,"_",All_clean$Domain)
All_clean$fasta<-paste0(">",All_clean$ID,"\n",All_clean$seq)


querydb<-All_clean$fasta[which(All_clean$hgnc_symbol!="CDH17")]
querydb<-All_clean$fasta
targetdb<-All_clean$fasta[which(All_clean$hgnc_symbol=="CDH17")]

write.table(data.frame(x=querydb),"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/mmseqs2/query.fasta",quote = F,col.names = F,row.names = F)
write.table(data.frame(x=targetdb),"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/mmseqs2/target.fasta",quote = F,col.names = F,row.names = F)


#> 2. Extract Extracellular FULL length sequence

CDHfamily_full<-data.frame()
for(jj in 1:length(ENTREZID$hgnc_symbol)){
  
  id<-ENTREZID$uniprotswissprot[jj]
  symbol<-ENTREZID$hgnc_symbol[jj]
  
  domaininfo<-GetProteinAnnontate(id,"ft_topo_dom")
  ori_seq<-GetProteinAnnontate(id,"sequence")
  
  if(domaininfo=="NA"){
    domainposi<-data.frame(symbol=symbol,start=c("NA"),end=c("NA"),
                           Location=c("NA"),Info=c("NA"),seq=c("NA"))
  }else{
    
    domain_list<-unlist(str_split(domaininfo,"; TO"))
    
    tmplist<-str_split(domain_list,";")
    
    tmplist<-data.frame(do.call("rbind", tmplist))
    if(dim(tmplist)[2]==2){
      tmplist<-data.frame(cbind(tmplist,data.frame(X3="")))
    }
    colnames(tmplist)<-c("Position","Location","info")
    domainposi<-data.frame(do.call("rbind",str_split(gsub("[A-Z]+_[A-Z]+ ","",tmplist$Position),"\\..")))
    domainposi<-cbind(domainposi,gsub(" \\/note=","",tmplist$Location))
    domainposi<-cbind(domainposi,tmplist$info)
    colnames(domainposi)<-c("start","end","Location","Info")
    domainposi$seq<-NA
    for(ii in 1:(dim(domainposi)[1])){
      domainposi$seq[ii]<-substr(ori_seq, start = as.numeric(domainposi$start[ii]), 
                                 stop = as.numeric(domainposi$end[ii]))
    }
    tmpc<-data.frame(symbol=rep(symbol,dim(domainposi)[1]))
    domainposi<-data.frame(cbind(tmpc,domainposi))
  }
  CDHfamily_full<-data.frame(rbind(CDHfamily_full,domainposi))
  print(jj)
}


All_full<-merge(ENTREZID,CDHfamily_full,by.x="hgnc_symbol",by.y="symbol",all.y=T)
write.xlsx(All_full,"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/CDHfamily_full_seq.xlsx")

All_clean<-All_full[-which(All_full$Location=="NA"),]
All_clean$Location<-gsub(" ","",All_clean$Location)
All_clean$ID<-paste0(All_clean$hgnc_symbol,"_",All_clean$start,"_",All_clean$end)
All_clean$fasta<-paste0(">",All_clean$ID,"\n",All_clean$seq)
write.xlsx(All_clean,"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/CDHfamily_full_seq_fasta.xlsx")


querydb<-All_clean$fasta[which(All_clean$Location=="Extracellular")]
targetdb<-All_clean$fasta[which(All_clean$hgnc_symbol=="CDH17" & All_clean$Location=="Extracellular")]

write.table(data.frame(x=querydb),"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/mmseqs2_full/query.fasta",quote = F,col.names = F,row.names = F)
write.table(data.frame(x=targetdb),"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/mmseqs2_full/target.fasta",quote = F,col.names = F,row.names = F)


# aln_full<-read.xlsx("/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/mmseqs2_full/BLAST/resultDB_s11.xlsx")
# all_tmp<-All_clean[,c("hgnc_symbol","uniprotswissprot","Location","ID","seq")]
# all_tmp<-all_tmp[which(all_tmp$Location=="Extracellular"),]
# colnames(aln_full)[1]<-"ID"
# 
# all_aln_full<-merge(aln_full,all_tmp,by="ID",all=T)
# 
# write.xlsx(all_aln_full,"/Users/qiuluting/Desktop/CTMbio/Project/2_Pipline/D2/CDH17_and_Family_domain/result/resultDB_aln_full_seq.xlsx")
















