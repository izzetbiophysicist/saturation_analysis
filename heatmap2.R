######################################################
### Reads AlphaMissense data and create a
### heatmap of the saturation mutagenesis
######################################################

setwd("~/artemis/figura_paper/")

# Install and load required packages
library(tidyverse)
library(readr)
library(colorspace)
library(plotly)

output <- "heat_map.png"
#input <- "~/artemis/structure/rhapsody_saturation/rhapsody-predictions.txt"
input <- "artemis_alphamiss.tsv"

pred.input <- read.table(input, sep='\t')
pred.input$V2 <- as.character(pred.input$V2)

####################################################################################################
# get sequence information
pred <- pred.input
tmp <- data.frame(V2=parse_number(pred$V2), V3=substr(pred$V2, 1, 1), V4=substr(pred$V2, nchar(pred$V2), nchar(pred$V2)), V6=pred$V3)
resnumbers <- unique(sort(tmp$V2))
seq.native <- NULL
for(i in resnumbers){
  native <- data.frame(V2=i, V4=unique(sort(tmp$V3[which(tmp$V2 == i)])),V6=0)
  tmp_seq_native <- as.character(native$V4)
  seq.native <- c(tmp_seq_native,seq.native)
}  
####################################################################################################
# generate data frame 
####################################################################################################
make_heatmap <- function(pred){
  to_plot <- data.frame(V2=NA, V4=NA, V6=NA)
  tmp <- data.frame(V2=parse_number(pred$V2), V3=substr(pred$V2, 1, 1), V4=substr(pred$V2, nchar(pred$V2), nchar(pred$V2)), V6=pred$V3)
  pred <- tmp
  sequence <- unique(sort(pred$V2))
  for(i in sequence){
    tmp_pred <- pred[which(pred$V2 == i),c(1,3,4)]
    native <- data.frame(V2=i, V4=unique(sort(pred$V3[which(pred$V2 == i)])),V6=0)
    tmp_pred <- rbind(tmp_pred, native)
    to_plot <- rbind(to_plot, tmp_pred)
  }
  to_plot$V6[which(is.nan(to_plot$V6))] <- 0
  to_plot <- na.omit(to_plot)
  names(to_plot) <- c("Resno","ResID","score")
  to_plot$ResID <- factor(to_plot$ResID,levels = sort(unique(to_plot$ResID),
                                                      decreasing = T))
out <- to_plot  
}
# run function and generate ggplot data frame
df <- make_heatmap(pred.input)

###################################################################################################3
# Plot function
###################################################################################################3
textcol <- "black" # text color
sep.scale <- 50 # residue X-axis labels separation 

# to obtain the correct separation in the x-axis
aux <- round(max(df$Resno)/sep.scale)
lab.x <-c(1,seq(sep.scale,aux*sep.scale-1, sep.scale))

p <- ggplot(df, aes(Resno, ResID, fill= score)) + 
  geom_tile(alpha = 0.65, color = 'gray7', linewidth = 0.1) +
  coord_fixed(ratio = 7) +  # Set the aspect ratio to 1 for squares
  xlab("Amino acid index") +
  ylab("Residue") +
  theme_grey(base_size=10)+
  theme(legend.position="right", legend.direction="vertical",
        legend.title=element_text(colour=textcol),
        legend.margin=margin(grid::unit(0, "cm")),
        legend.text=element_text(colour=textcol, size=7, face="bold"),
        legend.key.height=grid::unit(0.8, "cm"),
        legend.key.width=grid::unit(0.2, "cm"),
        axis.text.x=element_text(size=10, colour=textcol),
        axis.text.y=element_text(vjust=0.2, colour=textcol),
        axis.ticks=element_line(size=0.4),
        plot.background=element_blank(),
        panel.border=element_blank(),
        plot.margin=margin(0.7, 0.4, 0.1, 0.2, "cm"),
        plot.title=element_text(colour=textcol, hjust=0, size=14, face="bold"))+
  scale_x_continuous(expand = c(0, 0),breaks = lab.x, labels = lab.x) +
  scale_y_discrete(expand = c(0, 0))+
  scale_fill_continuous_diverging(palette = "Blue-Red 3", mid = 0.5,l1 = 30, l2 = 100, p1 = .9, p2 = 1.2)
p
ggsave(p, filename=output, height=6, width=12, units="in", dpi=200)

###################################################################################################3
# insert clinical data into heatmap
###################################################################################################3

mut.clinvar <- read.csv('clinvar_result.txt',sep='\t')
mut.clinvar <- mut.clinvar[which(mut.clinvar$Germline.classification %in% c('Pathogenic', 'Likely pathogenic')),]
clinical <- c()
for(entry in 1:nrow(mut.clinvar)){
  ## Pegando o primeiro elemento da lista de mutacoes, vem apenas a canonica
  clinical <- c(clinical, strsplit(mut.clinvar$Protein.change[entry],split=', ')[[1]][1])
}

clinical <- data.frame(V1=clinical)

read_clinical <- function(data){
  data$V1 <- as.character(data$V1)
  tmp <- data.frame(Resno=parse_number(data$V1), ResID=substr(data$V1, nchar(data$V1), nchar(data$V1)), score = rep(x = 0, length(data$V1)))
  tmp$ResID <- as.character(tmp$ResID)
  for(i in 1:length(tmp$Resno)){
    tmp$score[i] <- df$score[df$ResID == tmp$ResID[i] & df$Resno == tmp$Resno[i]]
  }
out <- tmp
}

am.mut.clinvar <- read_clinical(clinical)
###################################################################################################3
# plot com as mutações destacadas

# colorido pelo score
p <- p + geom_point(data = am.mut.clinvar, aes(fill=score), shape = 22, alpha = 1, color='gray7', size=2)
ggsave(p, filename=output, height=6, width=12, units="in", dpi=200)
# cor igual para localizar mais facilmente
p <- p + geom_point(data = am.mut.clinvar, fill='green', shape = 22, alpha = 1, color='gray7', size=2)
ggsave(p, filename=output, height=6, width=12, units="in", dpi=200)

ggsave(p, filename=output, height=6, width=12, units="in", dpi=200)


  
