library(tidyverse)
library(paletteer)
library(ComplexHeatmap)
library(circlize)

dat = read_tsv("2.Highfreq.PMD.pairwise.Aver.distance_ppm10corrected.tsv") %>%
  filter(PMD.direct.weighted!="-") %>%
  group_by(Prefix,PMD.abs) %>%
  summarise(Count=n()) %>%
  pivot_wider(names_from = Prefix,values_from = Count)

rank = read_tsv("2.Highfreq.PMD.pairwise.Aver.distance_ppm10corrected.tsv") %>%
  filter(PMD.direct.weighted!="-") %>%
  group_by(Prefix,PMD.abs) %>%
  summarise(Count=n()) %>%
  group_by(PMD.abs) %>%
  summarise(n=sum(Count)) %>%
  arrange(-n)

dat.plot = data.frame(row.names = dat$PMD.abs, dat[,-1])
dat.plot[is.na(dat.plot)] = 0
color = paletteer_c("grDevices::Heat",n = 10,direction = -1)

pdf("2.PMD.freq.pdf",wi=5,he=6.5)
Heatmap(dat.plot,
        show_row_names = T,
        show_column_names = T,
        row_names_gp = gpar(col="black"),
        row_names_max_width = unit(0.1,"mm"),
        row_names_side = "right",
        rect_gp = gpar(col = NA),
        row_title = "PMD",
        row_title_side = "right",
        col = color, na_col = "#ffffff",
        row_dend_width = unit(30, "mm"),
        heatmap_legend_param = list(title= "Count",
                                    title_position ="topleft", 
                                    legend_direction="horizon", 
                                    row_names_gp = gpar(fontsize = 8), 
                                    column_names_gp = gpar(fontsize = 10)),
        cluster_rows = T, 
        cluster_columns = T,
        cell_fun = function(j, i, x, y, width, height, fill) {
          if(dat.plot[i, j]!=0){
            grid.text(sprintf("%.0f", dat.plot[i, j]), x, y, gp = gpar(fontsize = 4))}
          
        }
)
dev.off()