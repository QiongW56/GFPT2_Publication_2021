# GFPT2_Publication_2021
This repository is for the figures created in the GFPT2 paper in 2021

## Figure 1C and Figure 1D
R version: 3.5.2; ComplexHeatmap version: 1.20.0 <br/>

### Figure 1C
Note: the genes in the heatmap also includes genes without significant pvalue less than 0.05. For each protein, the SILAC **median values** for one cell line were used for the plot.

### Figure 1D
Note: the Heatmap plotting R file was copied from the MetabolicTargetHeatmap and mordifed based on that; the genes in the heatmap also includes genes without significant pvalue less than 0.05 <br/>
1. Import LFQ data from "Proteomics_LFQ_summary_16.04.2018,xlsx" <br/>
2. Import SILAC data from "1 Proteomics_SILAC_summary_11.10.2018.xlsx" <br/>
3. import the metabolic markers from the paper "Markers.xlsx" <br/>
4. find the markers in both LFQ and SILAC <br/>
5. Log2 transform <br/>
6. z-score <br/>
7. Heatmap plotting <br/>
8. for details, refer to the Rfiles <br/>
9. For each protein, the SILAC **mean values** for one cell line were used for the plot. <br/>

## Figure 1E
R version: 4.0.0; dendextend_1.13.4; ggdendro_0.1-20; ggplot2_3.3.0 <br/>
Perseus version: 1.6.2.2

Note: How it was generated? <br/>
LFQ: <br/>
        1, find proteins in both LFQ and Literature using Protein ID in Excel at first then R later; <br/>
        2, log2(), replacing missing values and z-score in Perseus <br/>
        3, merge two datasets in Persues <br/>
        4, export from Persues, plot in R <br/>
(LFQ didn't filter anything and just replace the missing values directly, cause the normal distribution was ok) <br/>

C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\1 DundeeProteomicDataSet_14112017\iBAQ_TNBC\PerseusTextFile\PerseusOutputTable_2541 proteins in both iBAQ and literature_z score.txt -> 7 2541 Proteins in iBAQ and literature and with all detected_z score.xlsx <br/>

<br/>

SILAC: <br/>
         1, Use GeneName to match SILAC and Literature (SILAC has GeneName and ProteinName, but no ProteinID) <br/>
         2, find proteins in both SILAC and Literature using GeneName in Excel <br/>
         3, log2() in Perseus <br/>
         4, filter valid values, for Literature, 70% total valid values and for SILAC, 2/3 as valid <br/>
         5, replacing missing values and z-score in Perseus <br/>
         6, merge two datasets in Perseus <br/>
         7, export from Perseus, plot in R" <br/>

C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\1 DundeeProteomicDataSet_14112017\iBAQ_TNBC\iBAQ_SILAC_2\PerseusIntputTable_SILAC_Literature_TNBC.xlsx (you can find all the proteins with ENSG in the raw file) -> Literature.txt & SILAC.txt -> PerseusOutputTable_SILAC_Literature_TNBC.sps -> PerseusOutputTable_SILAC_Literature_TNBC_5533.txt -> SILAC_Literature_TNBC_5385.xlsx (delete all the NAs)


