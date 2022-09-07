
library(argparser, quietly=TRUE)
p <- arg_parser("")
p <- add_argument(p, "--input_variant_file", help="")
p <- add_argument(p, "--h19_geneanno_file", help="", default="")
p <- add_argument(p, "--CNN_result_directory", help="", default="")

argv <- parse_args(p)
input_variant_file = argv$input_variant_file
h19_geneanno_file = argv$h19_geneanno_file
CNN_result_directory = argv$CNN_result_directory

predicted_snp <- read.csv(input_variant_file,header=F)
predicted_snp <- as.character(predicted_snp$V1)

aa<-matrix(unlist(strsplit(predicted_snp,"_")), ncol = 4, byrow = T)
aa <- data.frame(aa)

geneanno <- read.csv(h19_geneanno_file)
aa$closest_distance_to_tss <- apply(aa, 1, function(x) min(abs(as.numeric(geneanno[geneanno$seqnames%in%x[1], ]$CAGE_representative_TSS)-as.numeric(x[2]))) )
aa <- aa[aa$closest_distance_to_tss < 20000, ]
aa$closest_distance_to_tss <- NULL
write.table(aa, file = paste0(CNN_result_directory,'predicted.vcf') ,row.names = F,col.names = F,sep = "\t",quote = F)

