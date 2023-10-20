source("DESeq2_help_function.R")

res = load_data_run_DESeq2("K562_1hr4sU_SAHA", "K562_1hr4sU_noSAHA", "total")
print(head(res))
write.csv(res, "totalRNA_60min4sU_SAHAvsno_DEseq2.csv")