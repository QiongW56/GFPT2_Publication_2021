# GFPT2_Publication_2021
This repository is for the figures created in the GFPT2 paper in 2021

## Figure 1A-B
Generated in .pptx files

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

## Figure 2
Note: <br/>
LFQ:
1. in Perseus ""VolcanoPlots.sps"", output tables: PerseusOutPutTable_X vs. Y
2. in R: ""VolcanoPlotForLFQ_X.R"", input the output tables from Perseus (around 2,700 proteins)
3. plot in R

SILAC:
1. Input output table from Perseus (in folder: .\MatchingLFQ&SILAC)
2. in R: ""VolcanoPlotForSILAC_X.R"", input the output tables from Perseus (around 5,100 proteins)
3. plot in R
4. SILAC: Median

Double check the labels on the volcano plot by comparing them to the “supplementary data 7 – significant proteins”, lots of genes are not in the dataset:
1.	It doesn’t have to be match, cause the “supplementary data 7” was reporting mean of LFQ and SILAC, the proteins have to be detected in both and in the same direction also being significant in both.
2.	LAMB3 in SILAC-EM, it has a p value a bit over 0.05, so not in the supplementary data
3.	PKP2 in SILAC-HE, not significant in LFQ (check on the raw data in “C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\MatchingLFQ&SILAC\5&6”)
4.	PRSS23 in LFQ-HM, not valid in SILAC
5.	POMC in SILAC-HM, not detected in LFQ
6.	KRT1 in SILAC-HM, not detected in LFQ
7.	DCD in SILAC-HM, not detected in LFQ
8.	PCSK1N in SILAC-HM, not detected in LFQ
9.	AKR1C1 in SILAC-HM, not detected in LFQ
10.	ALDH1A3 in SILAC-HM, not significant in LFQ

## Figure 3A-C
Treemap plots:
1. take genes from 'SupplementaryData7_DifferentlyExpressedProteins'
2. put genes for D492 vs. D492M and D492 vs. D492HER2 separately into Reactome (up- and down-regulated genes are together)
3. Used default setting
4. download results
5. work in excel (FDR < 0.05)
6. manually pick up the metabolic pathways
7. deleted the 'metabolism of protein' and 'metabolism of RNA' the big groups (take big space in treemap and give less information)
8. used Supplementary Data 6 but didnt delete proteins with SILAC p value > 0.05 (which = C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\IPA\ImportFile\Proteomics_LFQ&SILAC_X.xlsx files)

## Figure 3D
R version: 3.5.2; ComplexHeatmap version: 1.20.0 <br/>
Literature: Dihydropyrimidine Accumulation Is Required for the Epithelial-Mesenchymal Transition <br/>
note: this paper didn't post the raw data, so i could not know the expression of GFPT2 in each cell line (just mes and non-mes groups) <br/>

## Figure 3E-F
Note: <br/>
for the “profile” table, check on the file: <br/>
C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\SurvivalCurveGeneration_17.09.2019\Outputs\ Effict-Organized-For-Plotting.xlsx or Effect_ALL_Sign_Final_all_29Genes.xlsx <br/>
a.	all subtypes of breast patients <br/>
b.	all 29 genes for analysis <br/>
c.	20% of top or bottom patients for plotting <br/>

## Figure 4A-C
Gene expression (RNA or protein) bar plots <br/>

## Figure 4D-I
Boxplots - CDH1&CDH2 expression with siGFPT2 treatments <br/>
EMT markers: <br/>
1. used the results from the SQRDL_siRNA1 experiment <br/>
2. for D492M_CDH2, used the data from before (C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\4 Experiments documents_QIONG\14 Knockdowns Experiments\Results\8 KD_CDH1&2\KD_CDH1&2_Plotting") <br/>

## Figure 4J-L
Growth curves of D492, D492M, and D492HER2 with siGFPT2 treatments <br/>
Growth Curve (Proliferation): microscope --> pictures --> Macro API & ImageJ --> R <br/>

## Figure 4M
Migration rates with siGFPT2 treatments <br/>
Migration: output from IncuCyte <br/>

## Figure 4N-O
Invasion with siGFPT2 treatments <br/>
Invasion: Evos --> pictures --> ImageJ & Macro API --> R <br/>

## Figure 5
Note: <br/>
"CCLE: download GFPT2 data from CCLE website --> z-score and plotting in R <br/>
HMS: download data from HMS website --> organize in excel --> z-score and plotting in R <br/>
used Smith's paper to characterize the breast cell lines for both CCLE and HMS <br/>
Read ""README.txt"" in C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\ProteomicsManuscript_03.07.2019\NewFigures&Data_2020\1-MainFigures\Fig.10 <br/>
Patients: download data from TCGA --> Calculate in R --> plot in excel" <br/>

README.txt:
DPYD paper (with MMS signaures): 
It doesn't have expression data for genes in each cell line

TNBC-iBAQ paper: 
It doesn't have HER2 data, just TNBC and Luminal.
We can only use this for the dendrogram plotting because we can use the iBAQ data

So I combined HMS dataset with Smith's paper (cell line characterization) and CCLE dataset with Smith's paper:
HMS: less cell lines
CCLE: more cell lines
Smith's paper for classification

I can combine the TNBC-iBAQ and Smith's paper together for cell characterization, but they don't really agree with each other. 
But i replaced the "Basal B" with "Mesenchymal" for three cell lines "BT549","MDAMB231", "MDAMB157"








