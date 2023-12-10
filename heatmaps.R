######################################################
### Reads Rhapsody or AlphaMissense data and create a
### heatmap of the saturation mutagenesis
######################################################

# Install and load required packages
library(ggplot2)
library(readr)

output <- "~/artemis/structure/rhapsody_saturation/heat_map.png"
#input <- "~/artemis/structure/rhapsody_saturation/rhapsody-predictions.txt"
input <- "~/artemis/structure/artemis_alphamiss.tsv"

pred <- read.table(input, sep='\t')

make_heatmap <- function(type='rhapsody', pred, output){

if(type == 'rhapsody'){
to_plot <- data.frame(V2=NA, V4=NA, V6=NA)

sequence <- unique(sort(pred$V2))
for(i in sequence){
  tmp_pred <- pred[which(pred$V2 == i),c(2,4,6)]
  native <- data.frame(V2=i, V4=unique(sort(pred$V3[which(pred$V2 == i)])) ,V6=0)
  tmp_pred <- rbind(tmp_pred, native)
  to_plot <- rbind(to_plot, tmp_pred)  

}

to_plot$V6[which(is.nan(to_plot$V6))] <- 0

to_plot <- na.omit(to_plot)
}
  
if(type == 'alphamiss'){
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
}

# Heatmap 
png(file=output, width = 64, height = 48, units = "in", res = 300)
# Plot code
p <- ggplot(to_plot, aes(V2, V4, fill= V6)) + 
  geom_tile() +
  coord_fixed(ratio=1) +  # Set the aspect ratio to 1 for squares
  xlab("Amino acid index") +
  ylab("Residue") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

print(p)

# Save the plot
dev.off()

}


make_heatmap(type='alphamiss', pred, output)
