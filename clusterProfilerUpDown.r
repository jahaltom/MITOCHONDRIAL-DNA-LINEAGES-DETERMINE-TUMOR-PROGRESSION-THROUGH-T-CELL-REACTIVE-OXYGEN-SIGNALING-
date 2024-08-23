library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggridges)
#library(ggupset)
library(pathview)
#library(filesstrings)
library(pathview)

########GO#########

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
BiocManager::install(organism, character.only = TRUE,ask=FALSE)
library(organism, character.only = TRUE)




files=c("unique_NZB_TE_ACT_vs_NZB_TE_RES","unique_WT_TE_ACT_vs_WT_TE_RES",
        "NZB_TE_RES_vs_WT_TE_RES","NZB_TE_ACT_vs_WT_TE_ACT","NZB_TE_ACT_vs_NZB_TE_RES","unique_NZB_TE_RES_vs_WT_TE_RES","WT_TE_ACT_vs_WT_TE_RES")

for (file in files){

# reading in data from deseq2
df = read.csv(paste(file,"/",file,".DGE.txt" , sep=""), header=TRUE,sep="\t")

# we want the log2 fold change
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$Gene_name

# omit any NA values
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


    GOFunc <- function(GO)
    {

        gse <- gseGO(geneList=gene_list,
                     ont =GO,
                     keyType = "SYMBOL",
                     nPerm = 10000,
                     minGSSize = 3,
                     maxGSSize = 800,
                     pvalueCutoff = 0.05,
                     verbose = TRUE,
                     OrgDb = organism,
                     pAdjustMethod = "BH")

        require(DOSE)

        #Dotplot
        pdf(paste(file,"/",file,"_",GO,"._UP_Dotplot.pdf",sep=""),width=12,height=13)
        try(print(dotplot(gse, showCategory=gse$Description[gse$NES > 0][1:10], split=".sign",font.size = 18) + facet_grid(.~.sign)))
        dev.off()
        #Dotplot
        pdf(paste(file,"/",file,"_",GO,"._DOWN_Dotplot.pdf",sep=""),width=12,height=13)
        try(print(dotplot(gse, showCategory=gse$Description[gse$NES < 0][1:10], split=".sign",font.size = 18) + facet_grid(.~.sign)))
        dev.off()


        # #Encrichment Map
        pdf(paste(file,"/",file,"_",GO,"._UP_EncrichmentMap.pdf",sep=""),width=12,height=13)
        try( print(emapplot(pairwise_termsim(gse), showCategory =  gse$Description[gse$NES > 0][1:10],cex_label_category = 0.85)))
        # print(emapplot(pairwise_termsim(gse), showCategory = showCategory = 10,cex_label_category = 0.85))
        dev.off()
        # #Encrichment Map
        pdf(paste(file,"/",file,"_",GO,"._DOWN_EncrichmentMap.pdf",sep=""),width=12,height=13)
        try(print(emapplot(pairwise_termsim(gse), showCategory =  gse$Description[gse$NES < 0][1:10],cex_label_category = 0.85)))
        # print(emapplot(pairwise_termsim(gse), showCategory = showCategory = 10,cex_label_category = 0.85))
        dev.off()
        
        
        h=gse$Description[gse$NES < 0][1:10]
        t=gse$Description[gse$NES > 0][1:10]
        gse2=gse[gse$Description %in% c(h,t),]

        gse2$Description <- factor(gse2$Description, levels = gse2$Description)
        pdf(paste(file,"/",file,"_",GO,"._NES_Barplot.pdf",sep=""),width=12,height=13)
        try(print(ggplot(gse2, aes(x = Description, y = NES,fill=gse2$p.adjust)) +
                 geom_bar(stat = "identity")  +
                 coord_flip()   )     )            
        dev.off()       




    #Write the sig NES and Descriptions to csv.
        tempgse=gse[gse$p.adjust<=0.05,]
        tempgse=tempgse[,c("NES","Description")]
        write.table(tempgse,paste(file,"/",file,"_",GO,"._Sig_NES_Descriptions.csv",sep="") ,sep = ',',row.names = TRUE)




        # #Category Netplot
        # # categorySize can be either 'pvalue' or 'geneNum'
        # pdf(paste(file,"/",file,"_",GO,".Netplot.pdf",sep=""),width=20,height=25)
        # print(cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3,font.size = 7))
        # dev.off()

        # #Ridgeplot
        # pdf(paste(file,"/",file,"_",GO,".Ridgeplot.pdf",sep=""),width=15,height=20)
        # print(ridgeplot(gse) + labs(x = "enrichment distribution",font.size = 7))
        # dev.off()






    }
    try(GOFunc("BP"))
    try(GOFunc("CC"))
    try(GOFunc("MF"))



    ##########KEGG##########


# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
 # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
df = df[!duplicated(df[c("Gene_name")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$Gene_name %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)



    kegg_organism = "mmu"
    kk2 <- gseKEGG(geneList     = kegg_gene_list,
                   organism     = kegg_organism,
                   nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   keyType       = "ncbi-geneid")

if (nrow(kk2) !=0){
    #Dotplot
    pdf(paste(file,"/",file,"._UP_KEGG_Dotplot.pdf",sep=""),width=12,height=13)
    try(print(dotplot(kk2, showCategory = kk2$Description[kk2$NES > 0][1:10], title = "Enriched Pathways" , split=".sign",font.size = 18) + facet_grid(.~.sign)))
    dev.off()
    #Dotplot
    pdf(paste(file,"/",file,"._DOWN_KEGG_Dotplot.pdf",sep=""),width=12,height=13)
    try(print(dotplot(kk2, showCategory = kk2$Description[kk2$NES < 0][1:10], title = "Enriched Pathways" , split=".sign",font.size = 18) + facet_grid(.~.sign)))
    dev.off()


    #Encrichment Map
    pdf(paste(file,"/",file,"._UP_KEGG_EncrichmentMap.pdf",sep=""),width=12,height=13)
    try(print(emapplot(pairwise_termsim(kk2), showCategory = kk2$Description[kk2$NES > 0][1:10],cex_label_category = 0.85)))
    dev.off()
    #Encrichment Map
    pdf(paste(file,"/",file,"._DOWN_KEGG_EncrichmentMap.pdf",sep=""),width=12,height=13)
    try(print(emapplot(pairwise_termsim(kk2), showCategory = kk2$Description[kk2$NES < 0][1:10],cex_label_category = 0.85)))
    dev.off()
    
    
    
    h=kk2$Description[kk2$NES < 0][1:10]
    t=kk2$Description[kk2$NES > 0][1:10]
    kk22=kk2[kk2$Description %in% c(h,t),]
    kk22$Description <- factor(kk22$Description, levels = kk22$Description)
    pdf(paste(file,"/",file,"_KEGG_NES_Barplot.pdf",sep=""),width=12,height=13)
    try(print(ggplot(kk22, aes(x = Description, y = NES,fill=kk22$p.adjust)) +
             geom_bar(stat = "identity")  +
             coord_flip()  ))                    
    dev.off()     
    
    
    #Write the sig NES and Descriptions to csv.
    tempkk2=kk2[kk2$p.adjust<=0.05,]
    tempkk2=tempkk2[,c("NES","Description")]
    write.table(tempkk2,paste(file,"/",file,".KEGG_Sig_NES_Descriptions.csv",sep="") ,sep = ',',row.names = TRUE)

    


    # #Category Netplot
    # # categorySize can be either 'pvalue' or 'geneNum'
    # pdf(paste(file,"/",file,".KEGG_Netplot.pdf",sep=""),width=20,height=25)
    # print(cnetplot(kk2, categorySize="pvalue", foldChange=kegg_gene_list, showCategory = 3,font.size = 7))
    # dev.off()

    # #Ridgeplot
    # pdf(paste(file,"/",file,".KEGG_Ridgeplot.pdf",sep=""),width=15,height=20)
    # print(ridgeplot(kk2) + labs(x = "enrichment distribution",font.size = 7))
    # dev.off()

}

    #All KEGG Pathways
    #https://www.genome.jp/kegg/pathway.html

    # # Produce the 1st native KEGG plot (PNG)
    # dme <-pathview(gene.data=kegg_gene_list, pathway.id=kk2$ID[1], species = kegg_organism)


    # # #A bunch of mis KEGG plots
    # dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu01100 ",out.suffix="Metabolic pathways", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)
    # dme <- pathview(gene.data=kegg_gene_list, pathway.id="mmu00190 ",out.suffix="Oxidative phosphorylation", species = kegg_organism,kegg.native=FALSE,same.layer=FALSE)




 }



