#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("gridExtra")

library(dplyr) 
library(ggplot2)
library(reshape2)
library(gridExtra)


coverage_NC <- read.table("/home/lui4001/data/SAMPLE/output/aln_coverage.txt", header = FALSE, sep = "\t", col.names = c("chromosome", "depth", "numbases", "size", "percent_depth"))

genome <- coverage_NC %>% filter(grepl('genome', chromosome))

chr_index <- c("NC_000001", "NC_000002", "NC_000003", "NC_000004", "NC_000005", "NC_000006", "NC_000007", "NC_000008", "NC_000009", "NC_000010", "NC_000011", "NC_000012", "NC_000013", "NC_000014", "NC_000015", "NC_000016", "NC_000017", "NC_000018", "NC_000019", "NC_000020", "NC_000021", "NC_000022", "NC_000023", "NC_000024", "NC_012920", "genome")

id_index <- (c((seq(1:22)),"X", "Y","m", "total"))

avg_chr_cov <- vector()

for (i in seq(1:length(chr_index))){
  chro <- coverage_NC %>% filter(grepl(chr_index[i], chromosome))
  avg_chr_cov[i] <- sum(chro[,2]*chro[,3])/chro[1,4]
}

perc_chr_cov <- vector()
for (i in seq(1:length(chr_index))){
  chro <- coverage_NC %>% filter(grepl(chr_index[i], chromosome))
  perc_chr_cov[i] <- (1-chro[1,3]/chro[1,4])*100
}

summary <- data.frame(id_index, avg_chr_cov, perc_chr_cov)
summary$id_index <- factor(summary$id_index, levels = summary$id_index)

avg_cov_in <- ggplot(data=summary, aes(x=id_index, y=avg_chr_cov)) +
  geom_bar(stat="identity") + ggtitle("Average Coverage by Chromosome (w mtDNA)") +
  xlab("Chromosome ID") + ylab("Average Coverage")
per_cov_in <- ggplot(data=summary, aes(x=id_index, y=perc_chr_cov)) +
  geom_bar(stat="identity")+ggtitle("Percent Coverage by Chromosome (w mtDNA)") +
  xlab("Chromosome ID") + ylab("Percent Coverage")

avg_cov_ex <- ggplot(data=summary[c(1:24,26),], aes(x=id_index, y=avg_chr_cov)) +
  geom_bar(stat="identity")+ ggtitle("Average Coverage by Chromosome")+
  xlab("Chromosome ID") + ylab("Average Coverage")
per_cov_ex <- ggplot(data=summary[c(1:24,26),], aes(x=id_index, y=perc_chr_cov)) +
  geom_bar(stat="identity")+ggtitle("Percent Coverage by Chromosome") +
  xlab("Chromosome ID") + ylab("Percent Coverage")


summary_plots <- list(avg_cov_in, per_cov_in, avg_cov_ex, per_cov_ex)

avg_total <- summary$avg_chr_cov[26]
per_total <- summary$perc_chr_cov[26]

totals <- as.data.frame(cbind(avg_total, per_total))
colnames(totals) <-(c("Average", "Percent"))

write.csv(totals, "/home/lui4001/data/SAMPLE/output/totals.csv", row.names = FALSE)

ggsave("/home/lui4001/data/SAMPLE/output/coverage_summary.pdf", width = 11, height = 7.5, marrangeGrob(summary_plots, nrow=2, ncol=2))
