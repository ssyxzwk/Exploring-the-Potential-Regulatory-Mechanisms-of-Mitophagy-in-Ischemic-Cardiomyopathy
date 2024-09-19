## 2.1 Data Extraction


ICM-related datasets (GSE116250 and GSE5406) were acquired from the Gene Expression Omnibus (GEO) database (https://www.ncbi.nlm.nih.gov/geo/). The training set GSE116250 contained 13 ICM left ventricles samples and 14 controls (19). GSE5406 included 108 ischaemic cardiomyopathy samples and 16 controls and served as a validation set (20). The GSE116250 and GSE5406 datasets were chosen to investigate the potential regulatory mechanisms of mitochondrial autophagy in ischaemic cardiomyopathy based on the considerations that these datasets cover a wider patient population, provide more detailed clinical information, have a high level of data quality and reliability, and have a high degree of fit with the study aims and hypotheses. Selection criteria and steps for datasets: firstly, relevant datasets were screened according to the type of disease under study, ischaemic cardiomyopathy. Then, the sample sizes of the datasets were further evaluated, and preference was given to those datasets with sufficiently large sample sizes to support statistical analyses and testing of hypotheses. Finally, after comprehensively considering the rationality of experimental design, data availability, and research needs, GSE116250 and GSE5406 were selected as the final research subjects. Additionally, a total of 29 mitophagy related genes (MRGs) were mined from Pathway Unification database (https://pathcards.genecards.org/). A flowchart of the data analysis pipeline in this study was shown in Supplementary Figure S1.


## 2.2 Creation of Weighted Gene Co-expression Network
Using the ssGSEA algorithm of the “GSVA” R package (version 1.38.2) (21), the expression data of 29 MRGs were firstly extracted, and the enrichment scores of MRGs for each sample in the GSE116250 dataset were calculated. Subsequently, weighted gene co-expression network analysis was performed with the help of the “WGCNA” R package (version 1.70-3) (22), and the enrichment scores calculated above were used as sample features to construct and identify gene modules that were significantly co-expressed with MRGs. The samples were clustered to remove outliers. To ensure that gene interactions conform to the scale-free distribution to the maximum extent possible, we then performed soft threshold (β) determination. Topology selection for scale-free networks is more critical and has more stringent selection criteria compared to connectivity selection. Thereafter, the similarity among genes were calculated according to the adjacency, and the gene dendrogram including different modules was built using the degree of dissimilarity. Subsequently, the minimum number of genes per gene module was set to 100 and modules were merged using a correlation threshold of 0.2 to finalise the number of modules. Among these modules, we selected the module with the strongest correlation with MRGs enrichment scores (Cor = 0.44, P < 0.05) as the key module and obtained the final key module genes by plotting a heat map of the correlation between modules and traits.


## 2.3 Screening of Differentially Expressed Genes (DEGs) between ICM and control Group 


In order to identify genes with significant differences in expression between ICM and control groups and to further reveal the role of these genes in specific biological processes. The differentially expressed genes (DEGs) between ICM and control groups in GSE116250 were screened using “limma” R package (version 3.46.0) with |log2FC| > 1 and adj P < 0.05 (23). In addition to traditional microarray data analysis, limma has now been widely extended to RNA sequencing data analysis, where it excels in handling small samples and high-dimensional data, as well as controlling false positive rates. Then, we took the intersection between DEGs and MSRGs to gain differentially expressed MSRGs (DE-MSRGs) through ggVennDiagram (version 1.2.2). Subsequently, Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) were conformed using the “clusterProfiler” R package (version 4.4.4) (P < 0.05) and the "org.Hs.eg.db" R package (version 3.12.0) (P < 0.05) for up- and down-regulated DE-MSRGs (24). Besides, STRING database (https://cn.string-db.org/) was adopt to create the protein-protein interaction (PPI) network for these genes with confidence = 0.15.


## 2.4 Screening of the Biomarkers
The Least Absolute Shrinkage and Selection Operator (LASSO) (via “glmnet” R package (version 4.1-1)) and Support Vector Machine-Recursive Feature Elimination (SVM-RFE) were utilized to chose the characteristic genes, separately (25). After that, the intersection of two machine learning-obtained characteristic genes was established to gain key characteristic genes, which were defined as biomarkers. Afterwards, the diagnosis value of biomarkers for ICM was assessed and validated in GSE116250 and GSE5406 through Receiver Operating Characteristic (ROC) curves via “pROC” R package (version 1.17.0.1) (26). In addition, for the GSE5406 dataset, “PPROC” R package was employed to plot precision-recall (PRC) curves to further assess the performance of biomarkers in the diagnosis of ICM. The expression situation of the biomarkers between ICM and control groups were analyzed in GSE116250 and GSE5406, respectively. 
Additionally, based on the biomarkers, the nomogram was built to predict the risk rate of ICM via “regplot” R package (version 1.1) and “rms” R package (version 4.1-1), and corresponding calibration curve was drawn to verify the accuracy of nomogram.


## 2.5 Immune Analysis
In this study, we firstly calculated the content of 8 immune cells in all samples in GSE116250 using mcp counter algorithm of “immunedeconv” R package (version 2.0.4) (27). Then, the differences of content of immune cells between ICM and control groups were compared via Wilcoxon Test (P < 0.05). In order to understand the complex relationships and networks of interactions between immune cells, the relevant analysis of differential immune cells with each other was constructed (Cor > 0.3, P < 0.05). Also, spearman correlation analysis was performed for biomarkers and immune cells (Cor > 0.3, P < 0.05).

## 2.6 Construction of Regulatory Networks

To better explore the mechanisms involved in the biomarkers in ICM, we created the lncRNA-miRNA-mRNA network and transcription factor (TF)-mRNA network. The miRNAs regulating the biomarkers were predicted via miRTarBase database (https://mirtarbase.cuhk.edu.cn/), miRDB database (Http://www.mirdb.org), and TargetScan database (https://www.targetscan.org), and the predicted miRNAs of the three databases were intersected to obtain the common miRNAs. Then starbase database (Http://starbase.sysu.edu.cn) was utilized to predict the lncRNAs with regulatory interactions with common miRNAs (degraExpNum ≥ 1 and clipExpNum ≥ 1). Afterwards, The TFs targeting biomarkers were predicted by NetworkAnalyst database (https://www.networkanalyst.ca/). In the end, the Cytoscape software (version 3.8.2) was utilized to visualize the ceRNA network and TF-mRNA network (28).
2.7 Drug Prediction of the Biomarkers
Each biomarker was applied as key word to search drugs interacting with biomarkers in the DrugBank (https://go.drugbank.com/). The Cytoscape software (version 3.8.2) was used to visualize the biomarker-drug network (28).
2.8 Statistical Analysis
The R software (https://www.r-project.org/) was utilized to conduct statistical analysis. Differences analysis was executed via the Wilcox test (P < 0.05).


## lc.tableToNum ==  data[]<-lapply(data,as.numeric)