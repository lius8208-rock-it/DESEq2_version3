# DESeq2 version 3 is the third and last version of DESeq2 for differential gene expression analysis for my dataset. 
# The criteria of significantly differential expressed genes are adjusted p (padj)< 0.05 and log fold change (lfc) >2
# In the model, it only considered treatment and salinity, these two factors, no interaction term.
# The complied results for AW vs SW is actually only for coastal lake region, which is the served as baseline in this analysis. 
# comparing to version2, this version added likelihood ratio test (LRT) aiming to find consensus DEGs across three regions. The LRT was applied on interaction effect and regional effect. 
# More Venn diagram for the results comparing to Version2.