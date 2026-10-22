library(ggplot2)
library(ggpubr)
library(MASS)
library(dplyr)
library(tidyr)
library(broom)
library(coin)

data <- read.csv("OG.Branch_length.normalized.target.txt",sep = "\t",head = T)

ggplot() + geom_boxplot(data, mapping = aes(x = Node, y = norBL,fill=factor(Target)),outlier.shape=NA)+ylim(0,3)+
  scale_x_discrete(limits = c("N3", "N53", "N4", "N52", "N5", "N27", "N50", "N51", "N6", "N23", "N49", "N24", "N48", "N25", "N47", "N26", "N45", "N46", "N7", "N44", "N8", "N21", "N43", "N22", "N41", "N42", "N9", "N18", "N19", "N37", "N38", "N20", "N39", "N40", "N10", "N11", "N32", "N12", "N13", "N28", "N29", "N14", "N30", "N31", "N15", "N36", "N16", "N35", "N17", "N33", "N34"))
#scale_x_discrete(limits = c("N8","N21","N9","N18","N19","N20","N10","N11","N15"))

compare_MT_greater <- function(data, node_name) {
  sub_data <- data %>% filter(Node == node_name)
  mt <- sub_data %>% filter(Target == "QMT") %>% pull(norBL)
  nuc <- sub_data %>% filter(Target == "Nuclear") %>% pull(norBL)
  pt <- sub_data %>% filter(Target == "PT") %>% pull(norBL)
  results <- list()
  #MT vs Nuclear
  if (length(mt) >= 2 && length(nuc) >= 2) {
    test <- wilcox.test(mt, nuc, alternative = "greater")
    results[["MT_vs_Nuclear"]] <- tidy(test) %>% mutate(compare = "MT > Nuclear?")
  } else {
    results[["MT_vs_Nuclear"]] <- tibble(compare = "MT > Nuclear?", p.value = NA, 
                                         statistic = NA, estimate = NA, 
                                         note = "insufficient data")
  }
  #MT vs PT
  if (length(mt) >= 2 && length(pt) >= 2) {
    test <- wilcox.test(mt, pt, alternative = "greater")
    results[["MT_vs_PT"]] <- tidy(test) %>% mutate(compare = "MT > PT?")
  } else {
    results[["MT_vs_PT"]] <- tibble(compare = "MT > PT?", p.value = NA, 
                                    statistic = NA, estimate = NA,
                                    note = "insufficient data")
  }
  bind_rows(results) %>% mutate(Node = node_name, .before = 1)
}

all_nodes <- unique(data$Node)

all_results <- bind_rows(lapply(all_nodes, function(node) {
  compare_MT_greater(data, node)
}))
print(all_results)

significant_results <- all_results %>% filter(p.value < 0.005)
print(significant_results)

write.csv(all_results, "OG.Branch_length.normalized.target.wilcox.csv", row.names = FALSE)
