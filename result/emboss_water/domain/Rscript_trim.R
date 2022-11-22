setwd("/home/luting/emboss_water/domain")
query<-readLines("query.fasta")
#target<-readLines("target.fasta")
#ii<-1
output<-data.frame()

  for(ii in seq(1,length(query),2)){
    query_tmp<-query[c(ii,ii+1)]
    write.table(data.frame(x=query_tmp),"query_tmp.fasta",quote = F,col.names = F,row.names = F)
    system('sh cmd.sh',intern = TRUE)
    tmp<-read.table("result.tab",header=T,sep="\t", fill=T , row.names=NULL)
    colnames(tmp)<-names(tmp)[2:8]
    tmp<-tmp[,-8]
    # if(length(which(is.na(tmp$score))!=0)){
    #   for(nn in which(is.na(tmp$score))){
    #     tmp$start1[nn-1]<-tmp$row.names[nn]
    #     tmp$start2[nn-1]<-tmp$id2[nn]
    #   }
    #   tmp<-tmp[-which(is.na(tmp$score)),]
    # }
    # colnames(tmp)<-names(tmp)[2:12]
    # tmp<-tmp[,-12]
    output<-data.frame(rbind(output,tmp))
    print(query[ii])
}