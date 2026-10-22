library(ggridges)
library(ggplot2)

ks <- read.csv("all.ks.1022.tsv",sep="\t",header = T)

ks$Species <- factor(ks$Species, level = c("Thesium", "Santalum_Dendrotrophe", "Dendrotrophe", "Santalum", "Thesium_Santalum", "Erythropalum_Santalum", "Erythropalum", "Vitis", "Erythropalum_Vitis"))

ggplot(ks, aes(y = Species, x = Ks, fill = Species)) +
  geom_density_ridges_gradient() +
  theme_ridges() +
  labs(x = "Ks", y = "Species") +
  theme(legend.position = "None")+xlim(0,3)
