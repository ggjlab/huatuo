
#CNN_result_directory = '/home/ggj/QTL/Script_collating/Huatuo-NC/result/de_novo_predictions/'
#result_directory = '/home/ggj/QTL/Script_collating/Huatuo-NC/result/'
#Landscape_ieQTL_file = '/home/ggj/QTL/Script_collating/Huatuo-NC/data/Landscape_ieQTL.txt'
#coordinate_file='/home/ggj/hg19ToHg38.over.chain.gz'
#Landscape_ieQTL_Linked_SNPs_file = '/home/ggj/QTL/Script_collating/Huatuo-NC/data/Landscape_ieQTL_Linked_SNPs.RData'
#anno_file <- '/home/ggj/QTL/Script_collating/Huatuo-NC/data/HCL_Cell_Cluster_Annatation.txt'

library(argparser, quietly=TRUE)
library(dplyr)
library(data.table)
library(readr)

p <- arg_parser("")
p <- add_argument(p, "--CNN_result_directory", help="")
p <- add_argument(p, "--result_directory", help="", default="")
p <- add_argument(p, "--Landscape_ieQTL_file", help="", default="")
p <- add_argument(p, "--coordinate_file", help="", default="")
p <- add_argument(p, "--Landscape_ieQTL_Linked_SNPs_file", help="", default="")
p <- add_argument(p, "--anno_file", help="", default="")

argv <- parse_args(p)
CNN_result_directory = argv$CNN_result_directory
result_directory = argv$result_directory
Landscape_ieQTL_file = argv$Landscape_ieQTL_file
coordinate_file = argv$coordinate_file
Landscape_ieQTL_Linked_SNPs_file = argv$Landscape_ieQTL_Linked_SNPs_file
anno_file = argv$anno_file

#------------------------------------- Preprocess predicted results -------------------------------------#

file_list <- list.files(path = CNN_result_directory,pattern='no_gene.csv')
#file_list

cnn_result <- data.frame()
for (i in file_list){
   cnn_result_tmp <- fread(paste0(CNN_result_directory,i),sep=',',stringsAsFactors = F)
   cnn_result_tmp <- as.data.frame(cnn_result_tmp)
   cnn_result <- rbind(cnn_result,cnn_result_tmp)
  }
#cnn_result[1:3,1:10]

cnn_result$variant_ID <- paste0(cnn_result$'0','_',cnn_result$'1','_',cnn_result$'2','_',cnn_result$'3')
cnn_result <- cnn_result[!duplicated(cnn_result$variant_ID),]
cnn_result <- cnn_result[nchar(as.character(cnn_result$'2'))%in%1&nchar(as.character(cnn_result$'3'))%in%1, ]

cnn_result <- cnn_result[,c(7:ncol(cnn_result))]
#cnn_result[1:3,]

aa<-matrix(unlist(strsplit(cnn_result$variant_ID,"_")), ncol = 4, byrow = T)
variant_bed_hg19 <- data.frame(chr = aa[,1], start=as.numeric(aa[,2])-1, end = as.numeric(aa[,2]),
                              index = cnn_result$variant_ID)
options(scipen = 200)
write.table(variant_bed_hg19, file=paste0(CNN_result_directory,'variant_bed_hg19'), sep='\t', quote=F, row.names=F, col.names=F)

system(paste0("liftOver ", CNN_result_directory, "variant_bed_hg19 ", coordinate_file, " ",
              CNN_result_directory, "variant_bed_hg38 ", CNN_result_directory, "'NotMap'"),
       intern=FALSE, wait=TRUE)

variant_bed_hg38 <- read.table(paste0(CNN_result_directory,'variant_bed_hg38'), sep='\t')
bb <- matrix(unlist(strsplit(as.character(variant_bed_hg38$V4),"_")), ncol = 4, byrow = T)
variant_bed_hg38$variant_ID_h38 <- paste0(variant_bed_hg38$V1,'_',variant_bed_hg38$V3,'_',bb[,3],'_',bb[,4],'_b38' )
rownames(variant_bed_hg38) <- as.character(variant_bed_hg38$V4)

cnn_result <- cnn_result[cnn_result$variant_ID%in%rownames(variant_bed_hg38), ]
cnn_result$variant_ID_h38 <- variant_bed_hg38[as.character(cnn_result$variant_ID),]$variant_ID_h38
save(cnn_result , file=paste0(CNN_result_directory,'cnn_result.RData'))


#------------ Combine de novo prediction and variant-gene links from landscape-ieQTLs for each tissue -------------#

options(scipen = 5)
#load(paste0(CNN_result_directory,'cnn_result.RData'))
cnn_result_bei <- cnn_result

all_top_assoc_bei <- fread(Landscape_ieQTL_file,sep='\t',stringsAsFactors = F)
all_top_assoc_bei <- as.data.frame(all_top_assoc_bei)
all_top_assoc_bei$tissue <- gsub('\\d','',all_top_assoc_bei$celltype)
load(Landscape_ieQTL_Linked_SNPs_file)

tissues <- unique(gsub('\\d','',setdiff(colnames(cnn_result),c('variant_ID','variant_ID_h38'))))
tissues

combined_results <- data.frame()
for (tissue in tissues){
    
   #variant-gene links in a tissue
   all_top_assoc <- all_top_assoc_bei[all_top_assoc_bei$tissue%in%tissue, ]
   all_top_assoc$index <- paste0(all_top_assoc$variant_id , '-' , gsub('\\..*','',all_top_assoc$phenotype_id) )
   all_top_assoc <- all_top_assoc[order(all_top_assoc$index , all_top_assoc$pval_emt), ]
   all_top_assoc <- all_top_assoc[!duplicated(all_top_assoc$index), ]
   
   LD1 <- LD[LD$SNP_A%in%as.character(all_top_assoc$variant_id) , c('SNP_A','SNP_B','R2') ]
   colnames(LD1) <- c('variant_id','variant_ID_h38','R2')
   all_top_assoc <- merge(all_top_assoc,LD1,by="variant_id")
   
   #map ieqtl summary data of variant-gene links to prediction
   predict_model <- colnames(cnn_result_bei)[grep(colnames(cnn_result_bei), pattern=tissue)]
   cnn_result <- cnn_result_bei[cnn_result_bei$variant_ID_h38%in%as.character(all_top_assoc$variant_ID_h38), 
                                c('variant_ID_h38', predict_model) ]
   
   for (i in predict_model){
      tmp <- cnn_result[,c('variant_ID_h38',i)]
      tmp <- tmp[abs(tmp[,2]) > 0.5, ]
      if (nrow(tmp) > 0){
	     combined_results_tmp <- merge(all_top_assoc[,c('variant_id','phenotype_id','tss_distance', 'tests_emt',
                                                'pval_gi', 'pval_emt', 'pval_adj_bh','tissue','variant_ID_h38','R2')] 
                                                ,tmp,by="variant_ID_h38")
	     combined_results_tmp$predict_model <- i
         colnames(combined_results_tmp) <- gsub(i,'prediction',colnames(combined_results_tmp))
         combined_results <- rbind(combined_results , combined_results_tmp) 
	   }
    }   
}

colnames(combined_results) <- gsub('variant_id','top_ieqtl',colnames(combined_results))
save(combined_results , file = paste0(result_directory , 'combined_results.RData' ))


#---------------------------- cell-type-specific predicted variants -----------------------------#

cnn_result  <- cnn_result_bei
rownames(cnn_result) <- as.character(cnn_result$variant_ID_h38)
cnn_result <- cnn_result[,1:(ncol(cnn_result)-2)]
cnn_result <- cnn_result[rowSums(abs(cnn_result) > 0.5)>0,]

cnn_result_f <- reshape2::melt(as.matrix(cnn_result))
colnames(cnn_result_f) <- c('variant_ID_h38','predict_model','prediction')

#celltype_anno
celltype_anno <- read.csv(HCL_anno_file,sep='\t',header=F,stringsAsFactors=F) 
colnames(celltype_anno) <- c('Fined','Rough1')
celltype_anno$Fined1 <- gsub('.*\\(','',celltype_anno$Fined)
celltype_anno$Fined1 <- gsub('\\)','',celltype_anno$Fined1)
rownames(celltype_anno) <- celltype_anno$Fined1

cnn_result_f <- cnn_result_f[cnn_result_f$predict_model%in%rownames(celltype_anno), ]
cnn_result_f$celltype_anno_Rough <- celltype_anno[as.character(cnn_result_f$predict_model), ]$Rough1
cnn_result_f$variant_ID_h38 <- as.character(cnn_result_f$variant_ID_h38)
cnn_result_f_orig <- cnn_result_f

a <- reshape2::acast(cnn_result_f, variant_ID_h38~celltype_anno_Rough, value.var="prediction" , fun.aggregate = median, fill=NaN)
a[is.na(a)] <- 0
cnn_result_f <- reshape2::melt(a)
colnames(cnn_result_f) <- c('variant_ID_h38','celltype_anno_Rough','prediction')
cnn_result_f$celltype_anno_Rough <- as.character(cnn_result_f$celltype_anno_Rough) 
cnn_result_f$variant_ID_h38 <- as.character(cnn_result_f$variant_ID_h38) 
cnn_result_f_sig <- cnn_result_f[abs(cnn_result_f$prediction) > 0.5,  ] 

m <- reshape2::acast(cnn_result_f_orig, variant_ID_h38~predict_model, value.var="prediction")
m[is.na(m)] <- 0
m <- abs(m)

f1 <- function(x,i,celltype_anno=celltype_anno){
  cluster <- intersect(rownames(celltype_anno[celltype_anno$Rough1%in%i,  ]),colnames(m))
  if (length(cluster) > 1){
     t=t.test(x[cluster],x[setdiff(names(x),cluster)],alternative = "greater") 
	 return(as.numeric(t$p.value))
	 }
  else{
     t=t.test( x[setdiff(names(x),cluster)], mu=x[cluster], alternative = "less")
     return(as.numeric(t$p.value)) }
  }

f2 <- function(x,i,celltype_anno=celltype_anno){
  cluster <- intersect(rownames(celltype_anno[celltype_anno$Rough1%in%i,  ]),colnames(m))
  return( median(x[cluster]) > median(x[setdiff(names(x),cluster)]) )
 }
  
celltype <- unique(cnn_result_f_sig$celltype_anno_Rough)
celltype <- setdiff(celltype , 'unknown')

predicted_cell_specific_functional_variant_d <- list()
for (i in celltype){     
	 	
	 f1_celltype <- apply(m, 1, f1 , i=i,celltype_anno=celltype_anno)
     variant_1 <- names(f1_celltype[f1_celltype < quantile(f1_celltype,0.1)])
	 f2_celltype <- apply(m, 1, f2 , i=i,celltype_anno=celltype_anno)
	 variant_2 <- names(f2_celltype[f2_celltype==TRUE])
     specific_variant <- intersect(variant_1,variant_2)
 
	 tmp <- cnn_result_f_sig[cnn_result_f_sig$celltype_anno_Rough%in%i, ]
	 specific_variant <- intersect(specific_variant,tmp$variant_ID_h38 )
	 
	 tmp <- data.frame(celltype_anno_Rough=i,variant_ID_h38=specific_variant )
	 predicted_cell_specific_functional_variant_d <- rbind(predicted_cell_specific_functional_variant_d,tmp) 
 }
write.table(predicted_cell_specific_functional_variant_d, file = paste0(CNN_result_directory , 'predicted_cell_specific_functional_variant.txt'),
            sep='\t',quote=F,row.names=F,col.names=F)

predicted_cell_specific_functional_variant_d$index <- paste0(predicted_cell_specific_functional_variant_d$celltype_anno_Rough,'-',
                                                        predicted_cell_specific_functional_variant_d$variant_ID_h38)

combined_results <- combined_results[combined_results$predict_model%in%rownames(celltype_anno), ]
combined_results$celltype_anno_Rough <- celltype_anno[as.character(combined_results$predict_model), ]$Rough1
combined_results$predict_model <- celltype_anno[as.character(combined_results$predict_model), ]$Fined
combined_results$index <- paste0(combined_results$celltype_anno_Rough,'-',combined_results$variant_ID_h38)

cell_type_specific_gene_regulation <- combined_results[combined_results$index%in%as.character(predicted_cell_specific_functional_variant_d$index), ]

cell_type_specific_gene_regulation$index <- paste0(cell_type_specific_gene_regulation$variant_ID_h38,
                                            cell_type_specific_gene_regulation$celltype_anno_Rough,
                                            cell_type_specific_gene_regulation$phenotype_id,
                                            cell_type_specific_gene_regulation$tissue )
cell_type_specific_gene_regulation <- cell_type_specific_gene_regulation[order(cell_type_specific_gene_regulation$index,
                                                                 -cell_type_specific_gene_regulation$R2,
                                                                 cell_type_specific_gene_regulation$pval_adj_bh,
																 -abs(cell_type_specific_gene_regulation$prediction)
																 ), ]
cell_type_specific_gene_regulation <- cell_type_specific_gene_regulation[!duplicated(cell_type_specific_gene_regulation$index), ]


cell_type_specific_gene_regulation <- cell_type_specific_gene_regulation[,c('celltype_anno_Rough','variant_ID_h38','predict_model','prediction',
                                                                  'top_ieqtl','phenotype_id','pval_adj_bh')]

colnames(cell_type_specific_gene_regulation) <- c('Cell_type',
                                           'Putative_functional_regulatory_variant(cell-type-specific)',
										   'Prediction_model','Predicted_variant_effect',
										   'Linked_landscape_ieQTL(r2>0.8))',
										   'Putative_regulated_gene(ensemble_id)',
										   'pval_adj_bh(linked_landscape_ieQTL)')

write.table(cell_type_specific_gene_regulation, file = paste0(result_directory , 'cell_type_specific_gene_regulation.txt'),
            sep='\t',quote=F,row.names=F,col.names=F)








