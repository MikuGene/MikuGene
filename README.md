<h2>MikuGene: An integrated data analysis software package for bioinformatics beginners.</h2>
MikuGene Bioinformatics Ecological Community. ---- <i>Lianhao Song</i> (CodeNight)<br>
Newest Version: 1.2.1<br>
<br>
<i>@.@ What I've done:</i><br>
<br>
1. <b>Song L#</b>, Xie H#, Tong F, Yan B, Zhang S, Fu E, Jing Q and Wei L. Association of lnc-IL17RA-11 with increased radiation sensitivity and improved prognosis of HPV-positive HNSCC. <i>J Cell Biochem.</i> 2019;120:17438-17448.<br>
2. <b>Song L</b>, Zhang S, Duan C, Ma S, Hussain S, Wei L and Chu M. Genome-wide identification of lncRNAs as novel prognosis biomarkers of glioma. <i>J Cell Biochem.</i> 2019;120:19518-19528.<br>
3. <b>Song L#</b>, Zhang S#, Yu S, Ma F, Wang B, Zhang C, Sun J, Mao X and Wei L. Cellular heterogeneity landscape in laryngeal squamous cell carcinoma. <i>Int J Cancer.</i> 2020;147(10):2879-2890<br>
<br>
<i>@.@ What I'm interested in:</i><br>
<br>
1. Data analysis: ScRNA-seq; ATAC-seq; Bulk-seq; Spatial transcriptome ... <br>
2. Software production: On-line analysis platform; One-step data analysis software ... <br>
<br>
<i>@.@ How to contact me:</i><br>
<br>
Email 1: songlianhao233@gmail.com <br>
Email 2: 2743623823@qq.com <br>
<br>
<i>@.@ How to use MikuGene R package:</i><br>
<br>
<b><i>devtools::install_github("MikuGene/MikuGene")</i></b><br>
<br>
<i>@.@ On-line analysis platform: </i><br>
<br>
<b><i>http://hpvgroup.imwork.net/MKCell/</i></b><br>
<br><hr>

<h2>Functions:</h2>
<b><i>1. Single-cell RNA-seq data analysis:</i></b><br>
<b><i>1) MK_scRNA(x, name = "temp", Reso = 0.8, nGene = c(200,Inf), nVar = 2.5, Dim = 2, SCT = F, BatchRemove = F, Umap = F, Plot = T, Norm = T, save = T)</i></b><br>
 <b>x:</b> A two-dimensional matrix, including sparse matrices.<br>
 <b>name:</b> Custom naming of this process, such as file name for automatic backup storage, etc. Default 'temp'.<br>
 <b>Reso:</b> The resolution of cell clustering. Default 0.6. FindCluster() can be used to change the resolution of the output result.<br>
 <b>nGene:</b> Screen cells by setting an interval for the number of genes contained in each cell in standard procedure. Default c(200, Inf).<br>
 <b>nVar:</b> Filter the number of highly variable genes (k). Default 2.5 (2500 variable genes).<br>
 <b>Dim:</b> The dimension of unsupervised clustering. Default 2.<br>
 <b>SCT:</b> The 'sctransform' process recommended by Seurat will be used. Default 'SCT = F'.<br>
 <b>BatchRemove:</b> Harmony will be applied to remove batch effects according 'orig.ident'. Default 'BatchRemove = F'.<br>
 <b>Umap:</b> 'RunUMAP' function in Seurat will be applied. Default 'Umap = F'.<br>
 <b>Norm:</b> 'LogNormalize' will be applied in standard procedure. Default 'Norm = T'. (It does not affect SCT = T.)<br>
 <b>Plot:</b> Messages and pictures in the process will be displayed. Default 'Plot = T'<br>
 <b>Save:</b> The process will automatically create a backup folder under the working path, and store the resulting object as a name_backup.rds file. Default 'Save = T'<br>
 <b><i>OUTPUT:</i></b> This function will output a Seurat object whose cell clustering has been analyzed.<br>
 <br>
<b><i>2) MK_singler(x, ref = "HPCA", mode = "main", cluster = NULL, Cells = 10, name = NULL, Save = T)</i></b><br>
 <b>x:</b> A two-dimensional matrix, including sparse matrices.<br>
 <b>ref:</b> Ref-data, including HPCA (<i>HumanPrimaryCellAtlasData</i>), BPED (<i>BlueprintEncodeData</i>) and DICE (<i>DatabaseImmuneCellExpressionData</i>). Default "HPCA".<br>
 <b>mode:</b> Mode to identify cell, including main and fine. Default "main".<br>
 <b>cluster:</b> Cluster identities for each cell in <b>x</b>.<br>
 <b>Cells:</b> The number of cells (k) used to run separately. Default 10 (10,000 cells).<br>
 <b>name:</b> Custom naming of this process. Default 'temp'.<br>
 <b>Save:</b> The process will automatically create a backup folder under the working path, and store the resulting object as a name.csv file. Default 'Save = T'<br>
 <b><i>OUTPUT:</i></b> This function will output a pruned cell-label.<br>
 <br>
<b><i>3) MK_Spatial(Dir = getwd(), name = NULL, Reso = 0.6, Verbose = T, save = T)</i></b><br>
 <b>Dir:</b> waiting process ...<br>
 <br><hr>
<b><i>2. Gene biological functional enrichment analysis:</i></b><br>
<b><i>MK_Enrich(x, EnID = "temp", CutP = 0.01, Save = T, Wid = 8, Hig = 8.3)</i></b><br>
 <b>x:</b> genes, such as c('CD19', 'MS4A1', 'CD79A', ...).<br>
 <b>EnID:</b> Custom naming of this process, such as file name for automatic storage, etc. Default 'temp'.<br>
 <b>CutP:</b> P value cutoff.<br>
 <b>Save:</b> The process will automatically create a EnID folder, and store the resulting files. Default 'Save = T'.<br>
 <b>Wid:</b> The width of the figure when it is automatically saved. Default 8.<br>
 <b>Hig:</b> The height of the figure when it is automatically saved. Default 8.3.<br>
 <b><i>OUTPUT:</i></b> This function will output the gene enrichment results (GO, KEGG and Reactome).<br>
 <br><hr>
<b><i>3. WGCNA gene co-expression analysis:</i></b><br>
<b><i>MK_WG_Tom(x, name = "temp", nGene = 10000, Save = T)</i></b><br>
 <b>x:</b> A two-dimensional matrix, whose rows represent features (genes) and columns represent samples.<br>
 <b>name:</b> Custom naming of this process, such as file name for automatic backup storage, etc. Default 'temp'.<br>
 <b>nGene:</b> The number of genes will be used for co-expression network analysis.<br>
 <b>Save:</b> The process will automatically create a backup folder under the working path, and store the resulting object as a name_WGtom_backup.rds file. Meanwile the co-expression network diagram will be drawn and saved. Default 'Save = T'<br>
 <b><i>OUTPUT:</i></b> This function will output a list (WG_Tom) of five objects: <b>x</b>, MEs (column is each ME, row is each sample), Colors (color corresponding to each gene), Power (the power of WGCNA soft-threshold), TOM (the calculated topological matrix).<br>
 <br><hr>
<b><i>4. Save the large scRNA-seq matrix:</i></b><br>
<b><i>1) MK_toMMs(x, name = "temp", Cells = 10, verbose = F, HK_bm = F, Mito_rm = T, AC_rm = T, RP_rm = T, RPLS_rm = T, MIR_rm = T, ATP_rm = T, IGXV_rm = T)</i></b><br>
 <b>x:</b> A two-dimensional matrix, including sparse matrices.<br>
 <b>name:</b> Custom naming of this process, which creating a 'name' folder automatically. Default 'temp'.<br>
 <b>Cells:</b> The number of cells (k) used to store the matrix separately. Default 10 (10,000 cells).<br>
 <b>verbose:</b> Messages in the process will be displayed. Default 'verbose = F'.<br>
 <b>HK_bm:</b> Perform batch correction through housekeeping genes. Default 'HK_bm = F'.<br>
 <b>Mito_rm:</b> Delete mitochondrial genes. Default 'Mito_rm = T'.<br>
 <b>AC_rm:</b> Delete lncRNA (ACxxxx) genes. Default 'AC_rm = T'.<br>
 <b>RP_rm:</b> Delete pseudogenes genes. Default 'RP_rm = T'.<br>
 <b>RPLS_rm:</b> Delete ribosome-related genes. Default 'RPLS_rm = T'.<br>
 <b>MIR_rm:</b> Delete miRNA genes. Default 'MIR_rm = T'.<br>
 <b>ATP_rm:</b> Delete ATP-related genes. Default 'ATP_rm = T'.<br>
 <b>IGXV_rm:</b> Delete IGXV genes. Default 'IGXV_rm = T'.<br>
 <b><i>OUTPUT:</i></b> This function will automatically save the <b>x</b> in its own way.<br>
 <br>
<b><i>2) MK_reads(path, verbose = T)</i></b><br>
 <b>path:</b> The folder path <i>MK_toMMs</i> has output (When reading, please ensure that the folder has not been modified).<br>
 <b>verbose:</b> Messages in the process will be displayed. Default 'verbose = T'.<br>
 <b><i>OUTPUT:</i></b> This function will read the <i>MK_toMMs</i> output.<br>
 <br><hr>
<b><i>5. Virus Mapping:</i></b><br>
<b><i>1) MK_BuildVirusRef(version = "2020.3", OutVs = "default", verbose = T)</i></b><br>
 <b>version:</b> The data version of virusite.org. Default "2020.3".<br>
 <b>OutVs:</b> The SeqID of the virus genome not included in alignment. If none, input NULL. Default 'default'.<br>
 <b>verbose:</b> Messages in the process will be displayed. Default 'verbose = T'.<br>
 <b><i>OUTPUT:</i></b> This function will build virus-reference locally.<br>
 <br>
<b><i>2) MK_VirMap(path_r1, path_r2, name = "temp", maxMiss = 3, GTF = T, Pair = T)</i></b><br>
 <b>path_r1:</b> Fastq R1 (or .gz).<br>
 <b>path_r2:</b> Fastq R2 (or .gz).<br>
 <b>name:</b> Custom naming of this process. Default 'temp'.<br>
 <b>maxMiss:</b> Max miss-match in Mapping. Default 3.<br>
 <b>GTF:</b> Use GTF annotation files. Default 'GTF = T'.<br>
 <b>GTF:</b> If pair end. Default 'Pair = T'.<br>
 <b><i>OUTPUT:</i></b> This function will output the virus-mapping results.<br>
 <br><hr>
<b><i>6. Data sets:</i></b><br>
<i><b>1) Cellmarkers:</b> data('Cellmarker')</i><br>
<i><b>2) Cellstate (mainly from CancerSEA):</b> data('Cellstate')</i><br>
<i><b>2) Ligand and receptor (mainly from Ramilowski et al., 2015):</b> data('LigandReceptor')</i><br>
Ramilowski, J.A., T. Goldberg, J. Harshbarger, E. Kloppmann, M. Lizio, V.P. Satagopam, M. Itoh, H. Kawaji, P. Carninci, B. Rost, and A.R. Forrest. 2015. A draft network of ligand-receptor-mediated multicellular signalling in human. Nat Commun. 6:7866.
