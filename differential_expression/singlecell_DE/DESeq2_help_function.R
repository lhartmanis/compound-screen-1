library(anndata)
library(DESeq2)

load_data_run_DESeq2 = function(cond1, cond2, layer, conditionfile="conditionFileDNBPE.csv", filterfile="barcodes_with_sufficient_sequence_depth_20k.txt")
{
    ok_barcodes = scan(filterfile, what=character())
    conditiontable = read.csv(conditionfile)
    conditiontable = conditiontable[conditiontable$type == "Exon", ]
    conditiontable = conditiontable[conditiontable$condition %in% c(cond1, cond2), ]
    conditiontable = conditiontable[conditiontable$XC_DNBPE %in% ok_barcodes, ]   
    barcodes1 = as.vector(conditiontable[conditiontable$condition == cond1, "XC_DNBPE"])
    barcodes = as.vector(conditiontable[, "XC_DNBPE"])

    ad = read_h5ad("SAHA_DNB.annotatedData.h5ad")

    if(layer == "total") {
        df = ad$to_df(layer = "new") + ad$to_df(layer = "old")
    } else {
        df = ad$to_df(layer = layer)
    }
    df = df[rownames(df) %in% barcodes, ]
    condition = rownames(df) %in% barcodes1
    dds = DESeqDataSetFromMatrix(t(df), DataFrame(condition), ~ condition)
    dds = DESeq(dds, sfType="poscounts")
    res = results(dds)
    res
}