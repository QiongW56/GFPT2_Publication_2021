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
CCLE (https://portals.broadinstitute.org/ccle) <br/>
HMS: download data from HMS website --> organize in excel --> z-score and plotting in R <br/>
HMS LINCS (http://lincs.hms.harvard.edu/db/datasets/20348/) <br/>

used Smith's paper to characterize the breast cell lines for both CCLE and HMS <br/>
Molecular characterization of breast cancer cell lines through multiple omic approaches<br/>
SE Smith 2017<br/>
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5460504/<br/>
additional file 2<br/>

Read ""README.txt"" in C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\ProteomicsManuscript_03.07.2019\NewFigures&Data_2020\1-MainFigures\Fig.10 <br/>
Patients: download data from TCGA --> Calculate in R --> plot in excel" <br/>
https://www.cbioportal.org/study/summary?id=brca_metabric <br/>

README.txt: <br/>
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

## Figuer 6A
### Acidic analysis procedure
1. TargetLynx -> peak picking
2. Export "Complete Summary"
3. Put "temp" .txt file into "C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\Glycan precursor analysis\KDexperimentRepeat_16.07.2019\QIONG\TargetLynxFile"
4. Upload "temp" .txt file in Excel and named "temp.xlsx"
5. Open and run R file "AcidicNeg_KD_Met_Q_16.07.2019.R"
6. "AcidicNeg_Knockdown_Metabolomics_16.07.2019.xlsx" will be saved in "C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\Glycan precursor analysis\KDexperimentRepeat_16.07.2019\QIONG"
7. Replace all "first Scr sample" with the "second Scr sample" and saved in another .xlsx file "AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019"
8. Open and run R file "PlotingPercentageQ_CellType_17.07.2019.R" and "WT_Q_16.07.2019.bar.plot.R"
9. "AcidicNeg_Knockdown_Metabolomics_3Scr_template.xlsx" file used "AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019.xlsx" as a template to generate "Acetylaspartate
" plotting.

### Basic analysis procedure
1. TargetLynx -> peak picking
2. Export "Complete Summary"
3. Put "temp" .txt file into "C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\Glycan precursor analysis\KDexperimentBasic_03.12.2019\TargetLynx"
4. Upload "temp" .txt file in Excel and named "temp.xlsx"
5. Open and run R file "BasicNeg_Knockdown_Metabolomics_NOV2019.R" in Rfiles folder.
6. "BasicNeg_Knockdown_Metabolomics_NOV2019.xlsx" will be saved in "C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\Glycan precursor analysis\KDexperimentBasic_03.12.2019"
8. Open and run R file "KDexperiment_Basic_Plotting_NoScrButOtherKDs_10.12.2019.R" and other Rfiles to plot (need to recall more about this part)
9. Note: The R files for plotting were NOT dependent on "BasicNeg_Knockdown_Metabolomics_NOV2019.xlsx", they are related to files in 
C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\Glycan precursor analysis\KDexperimentBasic_03.12.2019\Metabolomics_KD_EMH_Basic_03.12.2019 <br/>

Note: <br/>
README in C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\ProteomicsManuscript_03.07.2019\NewFigures&Data_2020\1-MainFigures\Fig.7

## Figure 6B
UDP-GlcNAc levels in the three cell lines; basic mode data; Succinate IS <br/>

## Figure 6C
13C tracing experiments; normalized to the total metabolite level from basic mode (data used for figure 6B). <br/>

## Figure 6D-F
data from the acidic mode; <br/>
Internal standards: UDP-GlcNAc - phenylalanine IS; Glutamate - Glutamic acid IS; cystathionine - succinate IS <br/>

## Figure 7A
Metabolic profiling of cancer cells reveals genome-wide crosstalk between transcriptional regulators and metabolism<br/>
https://www.nature.com/articles/s41467-019-09695-9<br/>

Cellminer<br/>
https://discover.nci.nih.gov/cellminer/loadDownload.do<br/>

Metabolic profiling of cancer cells reveals genome-wide crosstalk between transcriptional regulators and metabolism<br/>
https://www.nature.com/articles/s41467-019-09695-9#additional-information<br/>

## Figure 7B
H<sub>2</sub>O<sub>2</sub> treatment on GFPT2 in the MDAMB231 cell line. <br/>

## Figure 7C-E
H<sub>2</sub>O<sub>2</sub> treatment on glutathione and GSH; GSH treatment on GFPT2 <br/>

## Figure 7F
Wide type expression of GSH in the three cell lines; acidic model; Glutamic Acid IS 

## Figure 7G-N
SQOR RNA expression - siGFPT2 treatment in four cell lines with two siRNAs.

## Figure 8A-D
Figures generated in excels. <br/>
Secretome analysis <br/>
EGF - GFPT2 <br/>
Insulin - GFPT2 <br/>
IGF1R higher in D492HER2 <br/>

## Figure 8E
IPA analysis of the phosphoproteomic dataset (D492 vs. D492HER2)

## Figure 8F
Perseus motif enrichment analysis for D492 vs. D492HER2; Perseus version 1.6.2.3

## Figure 8G-I
Figures generated in excels. <br/>
GSK3B_pS9 <br/>
GSK3B RNA level <br/>
GSK3B protein level <br/>

## Figure 8J-M
GSK3B knockdown - GFPT2 RNA expression

## Figure 9
Summary of the study - generated in excel.
        
## Supplementary Figure 1A
1.Used the “Peptide” excel sheet from the raw data and produced in Perseus <br/>
2.Generated in Perseus (version 1.6.2.2) and plotted in R <br/>

## Supplementary Figure 1B-C
Peptide length covered by LFQ and SILAC methods.

## Supplementary Figure 1D-F
1.	Genes identified and quantified in both LFQ and SILAC
2.	at least 2 out of 3 replicates in LFQ
3.	at least 2 out of 3 ratios in SILAC
4.	Log2FC were plotted
5.	no information about pvalue, so genes could be significant or not"

To generate the scatter plots
1. In Perseus&excel, find genes in both LFQ and SILAC, no matter for pvalues or directions because this is what you want to plot to see the correlation so you cannot select data manually
2. mean of LFQ and median of SILAC for the triplicates
3. import into R for plotting the scatter plot
4. marked the >=1 or <=-1 (deleted in the end)
5. labeled the top 50 genes (deleted in the end)

## Supplementary Figure 1G-H 
Note: <br/>
Heatmap: Median for SILAC not mean (you can check on the Perseus file “Annotations for E M HER2_SILAC_1.6.2.2.sps” in the folder “ProteinSignatureInCells”)
LFQ: 1841 proteins for plotting <br/>
SILAC: 2838 proteins for plotting

## Supplementary Figure 2
Heatmap of the EMT markers from the dbemt database along with the Reactome pathway analysis <br/>
EMT markers were downloaded from http://dbemt.bioinfo-minzhao.org/ <br/>
Genes in each cluster were put into Reactome (2020.4.18, Version 72)<br/>

Note: comparing to the figure 1C, there are only 12 out of the 26 genes in figure 1C are in this figure.

## Supplementary Figure 3A-C
Note: <br/>
Data used for GO annotation in Perseus v1.6.12.0: <br/>
C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/IPA/ImportFile – or “SupplementaryData6_DifferentlyExpressedProteins.xlsx"

1. load significant differently expressed genes in Perseus (only LFQ FDR < 0.05, no limit for SILAC p value)
2. use the whole proteome detected in SILAC as background
3. GO annotation
4. Fisher exact test, Benjamin-Hochberg FDR < 0.02
5. enrichment factor > 2 or 1.5 or 1.25 – just for easy-to-read

How the plot was generated:<br/>
EM&HE&HM_SILAC&LFQ_DataInput_07.02.2019.R <br/>
output: PerseusInPutTable_SILAC_EM&HE&HM_07.02.2019.txt <br/> 
EM&HE&HM_SILAC&LFQ_07.04.2019.sps in here (C:\Users\lenovo\OneDrive - Háskóli Íslands\PC-HI\5 ProteomicPaper\Figures&Tables in the paper\Annotations\Perseus\EM&HE&HM_06.02.2019_LFQ&SILAC_significant\Perseus_1.6.12.0) <br/> 
PerseusOutPutTable_X_Y_07.04.2019.txt <br/> 
EM&HE&HM_SILAC&LFQ_Plotting_07.02.2019.R for plotting




