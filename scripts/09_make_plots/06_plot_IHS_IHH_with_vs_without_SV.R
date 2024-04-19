# This script plots the following:
# Gene tracks
# iHS (original)
# iHH (original)
# iHS (SV removed and shifted)
# iHH (SV removed and shifted)
# For comparison reason, the bottom two plots are shifted back so that
# the positions are comparable
#
# library(rehh)
library(data.table)
library(tidyverse)

# Define inputs -----------------------------------------------------------
chrval <- "Superscaffold_chr8"
start_position <- 20580000
end_position <- 20750000
startSVposition <- 20678727
SVsize <- 8798

# Read inputs -------------------------------------------------------------
## Original -----
dt1a <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS.txt") %>% 
  filter(CHR == chrval & POSITION >= start_position & POSITION <= end_position)
dt1a$colour <- "#183059"
dt1a$colour[dt1a$LOGPVALUE > 6] <- "red"

dt1b <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/chr8_scanhh.txt") %>%
  filter(CHR == chrval & POSITION >= start_position & POSITION <= end_position)

dt1b <- dt1b %>% pivot_longer(cols = c("IHH_A", "IHH_D"), names_to = "IHH_type", values_to = "IHH")

## SV remove -----
dt2a <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/WGS_IHS_chr8_SVremoved_shifted.txt") %>%
  filter(CHR == chrval & POSITION >= start_position & POSITION <= (end_position - SVsize))
dt2a$POSITION[dt2a$POSITION > startSVposition] <- dt2a$POSITION[dt2a$POSITION > startSVposition] + 8798
dt2a$colour <- "#183059"
dt2a$colour[dt2a$LOGPVALUE > 6] <- "red"

dt2b <- fread("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/chr8_SVsnpsremoved_pos_shifted_scanhh.txt") %>%
  filter(CHR == chrval & POSITION >= start_position & POSITION <= (end_position - SVsize))
dt2b$POSITION[dt2b$POSITION > startSVposition] <- dt2b$POSITION[dt2b$POSITION > startSVposition] + 8798

dt2b <- dt2b %>% pivot_longer(cols = c("IHH_A", "IHH_D"), names_to = "IHH_type", values_to = "IHH")


ylab <- expression(paste(-log[10](italic("p")), " iHS", sep = ""))
p1a <- ggplot(data = dt1a) +
  geom_point(mapping = aes(x = POSITION/1000000, y = LOGPVALUE), colour = dt1a$colour, size = 0.5) + 
  scale_x_continuous(expand = c(0,0)) +
  ylab(ylab) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.y = element_text(size = 13))


p1b <- ggplot(dt1b) +
  geom_point(mapping = aes(x = POSITION/1000000, y = IHH, colour = IHH_type, pch = IHH_type)) + 
  scale_x_continuous(expand = c(0,0)) +
  ylab("iHH") +
  theme_bw() +
  scale_colour_manual(name = "iHH", values = c("red", "blue"), breaks = c("IHH_A", "IHH_D"), labels = c(expression("iHH"[Major]), expression("iHH"[Minor]))) +
  scale_shape_manual(name = "iHH", values = c(1, 2), breaks = c("IHH_A", "IHH_D"), labels = c(expression("iHH"[Major]), expression("iHH"[Minor]))) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.y = element_text(size = 13),
        legend.position = c(0.98,0.98),
        legend.background = element_rect(colour = "black", linewidth = 0.1),
        legend.justification = c(1,1)) 

p2a <- ggplot(data = dt2a) +
  geom_point(mapping = aes(x = POSITION/1000000, y = LOGPVALUE), colour = dt2a$colour, size = 0.5) + 
  scale_x_continuous(expand = c(0,0)) +
  ylab(ylab) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.y = element_text(size = 13)) 


p2b <- ggplot(dt2b) +
  geom_point(mapping = aes(x = POSITION/1000000, y = IHH, colour = IHH_type, pch = IHH_type)) + 
  scale_x_continuous(expand = c(0,0)) +
  ylab("iHH") +
  theme_bw() +
  scale_colour_manual(name = "iHH", values = c("red", "blue"), breaks = c("IHH_A", "IHH_D"), labels = c(expression("iHH"[Major]), expression("iHH"[Minor]))) +
  scale_shape_manual(name = "iHH", values = c(1, 2), breaks = c("IHH_A", "IHH_D"), labels = c(expression("iHH"[Major]), expression("iHH"[Minor]))) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.y = element_text(size = 13),
        legend.position = c(0.98,0.98),
        legend.background = element_rect(colour = "black", linewidth = 0.1),
        legend.justification = c(1,1)) 


# Create gene tracks for plotting -----------------------------------------
## Creat dummy gene track table for plotting ------------------------------
dt_genetracks <- data.table(chr = chrval,
                            POS = NA,
                            value = c(0.9,2.4),
                            data = "",
                            colour = NA)

## Read genes -------------------------------------------------------------
Genefile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/outlier_region/chr8_20500000_20800000_genes.txt"
dt_genes <- fread(Genefile, 
                  header = F)

dt_genes <- data.table(CHR = dt_genes$V1,
                       startPOS = dt_genes$V4,
                       endPOS = dt_genes$V5,
                       ID = gsub(pattern = "ID=", replacement = "", dt_genes$V9))

dt_genes <- dt_genes |> 
  filter(endPOS > start_position & endPOS < end_position)
dt_genes$starty <- c(1, 1, 1, 1, 1, 1, 1)
dt_genes$endy <- dt_genes$starty + 1
dt_genes$label <- 1:7
dt_genes$fill <- "chartreuse4"
dt_genes$fill[dt_genes$label == 6] <- "orange"
dt_genes$ylabpos <- 2

# Make plots ------------------------------------------------------------------
theme_plots <- theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y = element_text(face = "bold", colour = "black", size = 15),
        axis.title.x = element_text(face = "bold", colour = "black"),
        panel.grid = element_blank())

theme_empty <- theme_plots + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))
## Gene tracks -------------------------------------------------------------
x_axis_brks <- seq(20600000, end_position-50000, by = 50000)
px <- ggplot() + 
  geom_rect(dt_rep, mappin = aes(xmin=startPOS/1000000, xmax=endPOS/1000000,
                                 ymin= starty, ymax= endy), fill = dt_rep$fill, col = NA) +
  geom_rect(dt_genes, mappin = aes(xmin=startPOS/1000000, xmax=endPOS/1000000,
                                   ymin= starty, ymax= endy), fill = dt_genes$fill, col = NA) +
  geom_text(dt_genes, mapping = aes(x = (endPOS/1000000)+0.0015, y = ylabpos, label = label)) + 
  xlab("Position (Mb)") +
  scale_y_continuous(limits = c(0.9, 2.3), expand = c(0,0)) +
  scale_x_continuous(limits = c(start_position/1000000, end_position/1000000), expand = c(0,0), position = "top",  breaks = x_axis_brks/1000000) +
  theme_plots + 
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0, 5.5, 2.5, 5.5), "points"))


pcomb <- cowplot::plot_grid(px,
                            p1a,
                            p1b, 
                            p2a,
                            p2b,
                            align = "v", 
                            ncol = 1, 
                            # axis = "lr", 
                            rel_heights = c(1.7, 4,4,4, 4), labels = c("A)", "B)", "C)", "D)", "E)"))
# pcomb
png("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/IHS_IHH_chr8_with_vs_without_SV.png", width = 8.3, height = 8, units = "in", res = 600)
pcomb
dev.off()
