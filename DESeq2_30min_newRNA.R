source("/home/danielr/SAHA_and_compound_screen/src/DESeq2_help_function.R")

res = load_data_run_DESeq2("K562_30min4sU_SAHA", "K562_30min4sU_noSAHA", "new")
print(head(res))
write.csv(res, "/home/danielr/SAHA_and_compound_screen/SAHA_singlecell2/newRNA_30min4sU_SAHAvsno_DEseq2.csv")