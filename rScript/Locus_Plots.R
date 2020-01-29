###
#
# Tn_Seq pipeline 4 : Make some graphs
#
# Author : Fx_stubbe
#
###

# 0) Libraries
# -------------------------------- #

library(tidyverse)
library(data.table)
library(ggthemes)
library(patchwork)

# 1) Paths and Inputs
# -------------------------------- #

Data <- read_delim("/Users/carlo/Desktop/Eme/Output/full_pos.txt", delim = "\t")
gff <- read_csv("/Users/carlo/Desktop/Eme/Data/Melitensis16M_gene_table_v3.csv")


# 2) Fucntions
# -------------------------------- #

splitWithOverlap <- function(vec, seg.length, overlap) {
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}


# 2) Get a gene
# -------------------------------- #

#Old Locus NAme
Gene <- "BME_RS00795"
# Graph Boundaries
boundaries <- 200
#Plot_by
S_factor <- 10

#Get the desired gene
idx <- which(gff$GenBank_locus_tag == Gene)

#Extract Position
g_start <- gff$start[idx] - boundaries
g_end <- gff$end[idx] + boundaries

#Make DF
Data_gene <- Data %>% filter(Pos %in% g_start:g_end) 
Data_gene <- splitWithOverlap(Data_gene$log, S_factor, 0) 
Data_gene <- rapply(Data_gene, mean)
Data_gene <- tibble(Pos = 1:length(Data_gene), log = Data_gene)

#Get standard deviation
Data_gene$minlog <- Data_gene$log - sd(Data_gene$log, na.rm = T)
Data_gene$maxlog <- Data_gene$log + sd(Data_gene$log, na.rm = T)

#---------Make plot
#Coverage
xtick<-seq(0, nrow(Data_gene),50)

p1 <- ggplot(Data_gene, aes(x = Pos, y = log)) +
  geom_line(color = "firebrick", size = 0.7) +
  labs(x = "", y = "Log10", caption = paste("Genomic locus + ",boundaries," nt neighborhood",sep = "")) +
  ylim(c(0, 5))  +
  geom_vline(xintercept = 20, linetype="dashed", 
               color = "grey40", size=0.8) + 
  geom_vline(xintercept = nrow(Data_gene) - 20, linetype="dashed", 
             color = "grey40", size=0.8) +
  scale_x_continuous(breaks=xtick) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.caption = element_text(hjust = 0))


#Coverage + standard deviation
p2 <- ggplot(Data_gene, aes(x = Pos, y = log)) +
  geom_ribbon(aes(ymin = minlog, ymax = maxlog), alpha = 0.5,
      fill = "darkseagreen3", color = "transparent") +
  geom_line(color = "aquamarine4", lwd = 0.7) +
  labs(x = "", y = "") +
  ylim(c(-2, 6))  +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
            
#Plot + smoothing            
p3 <- ggplot(Data_gene, aes(x = Pos, y = log)) +
  geom_point(color = "gray40", alpha = 0.3) +
  labs(x = "", y = " ", caption = "B. melitensis bv. 1 str. 16M : GCA_000007125.1 ASM712v1") +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 1000),
  #            se = F, size = 1.3, aes(col = "1000")) +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 100),
  #            se = F, size = 1, aes(col = "100")) +
  stat_smooth(method = "gam", formula = y ~ s(x, k = 10),
     se = F, size = 0.8, aes(col = "10")) +
  scale_color_manual(name = "k", values = c("darkorange2", "firebrick","dodgerblue3"))+
  theme_bw() + 
  ylim(c(0, 5))  +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.caption = element_text(),
        legend.position = 'none')
            
#---------Wrap te graph

p_full <- p1+p2 / p3
p_full <- wrap_elements(p_full) + ggtitle(paste(Gene," | ",gff$Product[idx], " | Orientation ", gff$strand[idx] , sep=""))  

pdf(paste("/Users/carlo/Desktop/Eme/Output/Graphs/", Gene ,".pdf", sep = ""), width = 15, height = 4)
print(p_full)
dev.off()

