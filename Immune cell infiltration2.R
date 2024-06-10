rm(list = ls())



data_gene=read.table(" ",sep="\t",header = T,check.names = F,row.names = 1)



data=data_gene[c(" "," "," "," "),]





immune=read.table("CIBERSORT-Results.txt",sep="\t",header = T,check.names = F,row.names = 1)
rownames(immune)
colnames(data)



data=data[,match(substr(rownames(immune),1,10),colnames(data))]
identical(substr(rownames(immune),1,10),colnames(data))


scat=cbind(immune,as.data.frame(t(data)))

cell=as.character(colnames(scat)[1:19])
pathway=as.character(colnames(scat)[20:23])

tex <- as.data.frame(matrix(0,ncol=4,nrow=19))
colnames(tex) <- c("pathway","celltype","r", "p")

path=c("0","0","0","0")
for (j in pathway) {
  
  x=scat[,j]
  tex$celltype <- cell
  for (i in cell) {
    
    y=scat[,i]
    
    
    r= format(cor(x,y),digits = 4)
    p= format(cor.test(x,y)$p.value,digits = 4)
    
    tex[match(i,tex$celltype),c(1,3,4)]=c(j,r,p)
    
  }
  
  path=rbind(path,tex)
  
}

d<-as.data.frame(path)
d<-d[-1,]


path=d
path$r=as.numeric(path$r)
path$p=as.numeric(path$p)
colnames(path)=c("gene","celltype","correlation","pvalue")




path$stars=ifelse(path$pvalue< 0.001,"***",ifelse(0.001< path$pvalue & path$pvalue < 0.01,"**",ifelse(0.01 < path$pvalue & path$pvalue < 0.05,"*"," ")))



library(ggplot2)

ggplot(path, aes(gene,celltype)) +
  geom_tile(aes(fill=correlation)) +
  geom_text(aes(label=stars), color="black", size=4) + 
  scale_fill_gradient2(low='blue', high='red',mid = 'white', limit=c(-1,1),  name=paste0("*    p < 0.05","\n\n","**  p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation")) + # 把P值添加到图例，这一步非常巧妙
  labs(x=NULL,y=NULL) + 
  theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
        axis.text.y = element_text(size=8,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank()) 


head(path)
library(ggplot2)

ggplot(path,aes(x=celltype,y=gene))+
  geom_point(aes(size=-log10(pvalue),
                 color=correlation))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_colour_gradient2(low = "blue",mid = "white", high ="red")+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3))






