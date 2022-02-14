library(UpSetR)
movies <- read.csv(system.file("extdata", "movies.csv", package = "UpSetR"), 
                   header = T, sep = ";")


upset(movies, nsets = 6, number.angles = 30, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Genre Intersections", sets.x.label = "Movies Per Genre", 
      text.scale = c(1.3, 1.3, 1, 1, 2, 0.75))


#converting dataframe of DE genes into format for upsetR
within=read.table("within-family-analysis_corrected.csv",sep=",",header=TRUE,stringsAsFactors=T)
#drop column29
#within = subset(within, select=-c(N029))

#ended up keeing 29 easier to just drop genes where totla is less than 1 - then its unique to that family dont care 
save=do.call(rbind,lapply(within, function(x) table(factor(x, levels=c(levels(unlist(within)))))))
#quick manual edits remoe col V1 
#remove the column with the empty counts
#name the first column with Names 
write.csv(save, file="withinfamily_formatted.csv")


within=read.csv("withinfamily_formatted.csv",sep=",",header=TRUE,stringsAsFactors=T)

dat= as.data.frame(x = t(within), stringsAsFactors = FALSE)
names(dat) <- dat[1,]
dat <- dat[-1,]

dat=as.data.frame(dat)
df2 <- mutate_all(dat, function(x) as.numeric(as.character(x)))

upset(df2,nsets=13,mb.ratio = c(0.40, 0.60),    
      text.scale = c(1.3, 1.3, 1, 1, 2, 2))
      




blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

df <- data.frame(
  group = c("TSC1", "TSC2", "both"),
  value = c(10, 17, 5)
)

bp<- ggplot(df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)

library(scales)
#piechart for TSC mutations 
pie + scale_fill_brewer("Blues") + blank_theme +
  theme(axis.text.x=element_blank())
  


library(ggplot2)
setwd("C:/Users/Owner/OneDrive - Arizona State University/Documents/TSCdata")
bookings=read.csv("piedata.csv")

data=as.data.frame(table(bookings$all))

data$percent=paste(round((data$Freq/sum(data$Freq)), digits=1)*100), "%")


data$percent=round((data$Freq/sum(data$Freq)),digits=4)

label = factor(data$Var1)

data$label <- factor(paste(data$Var1, data$percent), levels = paste(data$Var1, data$percent)) 


ggplot(data, aes(x = "", y = Freq, fill=Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()) +
  
  theme(panel.background = element_blank()) 
###only use this part when you want to see the percetnages and then apply them yourself 

+

 geom_label(aes(label = percent),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE)













#dont need
###############################################################################
#cut out genes where totla is less than 1 - then its unique to that family dont care 
within <- subset( within, select = -Names )
dd1 = within[,colSums(within) > 1]
#rename ENSG with gene names 
library("AnnotationDbi")
library("org.Hs.eg.db")
#RUN THIS FOR EACH RESSIG COMBINATION - will work on getting rid of the NAs later 
genes=colnames(dd1)
colnames(dd1) <- mapIds(org.Hs.eg.db,
                        keys=genes,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

#drop columns that dont have a name
dd1=dd1[!is.na(names(dd1))]
m = make_comb_mat(dd1)
upset(m)
#nsets = number of columns in dd1
upset(dd1,nsets=36,mb.ratio = c(0.20, 0.80),show.numbers = "no")
upset(dd1,nsets=50,mb.ratio = c(0.30, 0.70))
#run this to see which gene names are which ENSG - then use within-family-analysis_edit_dropunaffected.csv
#to controlF which insersects in the barplot 
mapIds(org.Hs.eg.db,
       keys=genes,
       column="SYMBOL",
       keytype="ENSEMBL",
       multiVals="first")





#dont use 

ggplot(data, aes(x = "", y = Freq, fill=Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank()) +
  scale_fill_manual(values = c("dodgerblue4","dodgerblue3","dodgerblue2","dodgerblue1","deepskyblue4",
                               "deepskyblue3","deepskyblue2","deepskyblue1",
                               "cadetblue4","cadetblue3","cadetblue2","cadetblue1","gray45")) +
  theme(panel.background = element_blank()) 