# Title     : FDR
# Objective : FDR
# Created by: ldy93
# Created on: 2020-10-06

table <- read.csv('./PositiveSelection/Results/summary/positive_selection_result.csv', header=T)
b <- (table$p_value)
c <- p.adjust(b, 'fdr')
table$q_value <- c
write.csv(table, file='./PositiveSelection/Results/summary/positive_selection_adjusted.csv')
