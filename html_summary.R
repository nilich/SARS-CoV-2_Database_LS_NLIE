##### Script to test Table display for SARS2 Pipeline
library(kableExtra)
library(formattable)
library(dplyr)
library(knitr)

t <- read.csv("/mnt/nfs/bio/Sequencing/Virology/SARS-CoV-2_Database_LS_NLIE/Results/Summary_Results.csv", sep = ",", header=T)

t$Vadr_QC <- as.character(t$Vadr_QC)

t$X.Breath_10<- as.numeric(t$X.Breath_10)
t$X.Breath_400<- as.numeric(t$X.Breath_400)
row.names(t) <- NULL

t$Vadr_QC <- ifelse(
  t$Vadr_QC == "PASSED",
  cell_spec(t$Vadr_QC, color = "black"),
  cell_spec(t$Vadr_QC, color = "red")
)

t$X.Breath_10 <- ifelse(
  t$X.Breath_10 > 95,
  cell_spec(t$X.Breath_10, color = "black"),
  cell_spec(t$X.Breath_10, color = "red")
)

t$X.Breath_400 <- ifelse(
  t$X.Breath_400 > 90,
  cell_spec(t$X.Breath_400, color = "black"),
  cell_spec(t$X.Breath_400, color = "red")
)

t$MappedRads <- ifelse(
  t$MappedRads > 800000,
  cell_spec(t$MappedRads, color = "black"),
  cell_spec(t$MappedRads, color = "red")
)

colnames(t) <- c("Name", "Linage", "Propability", "MeanCoverage", "minCoV10x", "%minCoV10x", "minCoV400x", "%minCoV400x", "MappedReads", "TotalReads", "Vadr_QC")

kbl(t, escape = F, align = "l")%>%
  kable_styling(bootstrap_options = c("striped", "hover", "responsive"), full_width = F, position = "center") %>%
  save_kable(file = "test.html", self_contained = T)
