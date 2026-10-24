#PRUEBA CON UN VCF
#Anotado con Mane_select
library(optparse)
library(dplyr)
library(tidyr)  
library(GenomicRanges)

#Archivos VCF : Mejor poner direcciones de VCF en un txt o la direccion de la carpeta con los txt
option_list<-list(
  make_option("--vcf",help="text file containing list of vcf files, or directory containing all vcf files",dest='vcf'),
  make_option("--MANE",help='MANE_Select file to be used',dest='MANE'),
  make_option("--genes",help='gene panel file (first column = gene symbol, # = comment)',dest='genes'),
  make_option("--exclude",default="",help='comma-separated genes removed from the panel',dest='exclude')
)
opt<-parse_args(OptionParser(option_list=option_list))

vcf_file=opt$vcf 
MANE_file=opt$MANE
#
# Verificar si el archivo/directorio existe
if (!file.exists(vcf_file)) {
  stop("El archivo o directorio especificado no existe: ", vcf_file)
}

# Verificar si es un directorio de forma más robusta
if (dir.exists(vcf_file)) {                                                          
  filenames <- list.files(vcf_file, pattern = ".vcf", full.names = TRUE)            
  names <- list.files(vcf_file, include.dirs = FALSE, pattern = "+segments.vcf")
} else if (file.exists(vcf_file)) {                                                  
  filenames <- apply(read.table(vcf_file), 1, toString)                            
} else {
  stop("No se puede determinar si es archivo o directorio: ", vcf_file)
}

#Abrir el archivo MANE_select y armar el GRanges
MANE_file <- read.delim(MANE_file, sep = "\t", 
                        header = TRUE)
MANE_G <- makeGRangesFromDataFrame(MANE_file, 
                                 seqnames.field = "chr", 
                                 start.field = "start", 
                                 end.field = "end",
                                 keep.extra.columns = TRUE)

#Lista de genes: primera columna del panel germinal, menos los excluidos
if (is.null(opt$genes) || !file.exists(opt$genes)) {
  stop("Falta el panel de genes (--genes): ", opt$genes)
}
panel_lines <- trimws(readLines(opt$genes))
panel_lines <- panel_lines[panel_lines != "" & !startsWith(panel_lines, "#")]
genes_list <- unique(trimws(sub("\t.*", "", panel_lines)))
exclude <- trimws(strsplit(opt$exclude, ",")[[1]])
genes_list <- setdiff(genes_list, exclude)
if (length(genes_list) == 0) {
  stop("La lista de genes quedó vacía: ", opt$genes)
}
missing_mane <- setdiff(genes_list, MANE_file$gene_symbol)
if (length(missing_mane) > 0) {
  stop("Genes del panel ausentes en MANE: ", paste(missing_mane, collapse = ", "))
}

#Si el vcf_file es un directorio, se leen todos los archivos vcf dentro de ese directorio

if(file.info(vcf_file)$isdir){                                                          #if the input vcf file is a directory;
  filnames<-list.files(vcf_file,pattern=".vcf",full.names=T)                            #puts all .vcf files in that directory on the 'vcf' list
  names <- list.files(vcf_file,  include.dirs = FALSE, pattern = "+segments.vcf")
}else{                                                                                  #else expects vcf file to be a file containing a list of all the vcf files to be used
  filenames<-apply(read.table(paste(vcf_file)),1,toString)                            #and reads in that list
}
output_file <- filenames

#Se crea una lista de los vcfs que se van a 
vcf <- list()
vcf_G <- list()
for (i in 1:length(filenames)) {
#Lectura del archivo vcf  
  vcf[[i]] <- read.delim(filenames[i],  sep = "\t", 
                         header = FALSE, 
                         comment.char = "#",
                         stringsAsFactors = FALSE)
#Se le asignan los nombres de las columnas al vcf 
  names[i] <- sub(".vcf.*", "", names[i])
  colnames(vcf[[i]]) <- c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT",	names[i])
#Se extrae la columna INFO y se crea una nueva columna END con el valor de END del vcf
    vcf[[i]]$END <- sub(".*END=", "", vcf[[i]]$INFO)
  vcf[[i]]$END <- as.integer(vcf[[i]]$END)
#Se arman columnas disntintas para los paramentros del vcf  
  vcf[[i]] <- vcf[[i]] %>%
    separate(names[i], 
             into = c("GT", "CN", "NP", "QA", "QS", "QSE", "QSS"), 
             sep = ":")
#Convertir el data frame en Granges    
  vcf_G[[i]] <- makeGRangesFromDataFrame(vcf[[i]], 
                                         seqnames.field = "CHROM", 
                                         start.field = "POS", 
                                         end.field = "END",
                                         keep.extra.columns = TRUE)
#Superponer con MANE_select para anotar los exones
  res <- findOverlaps(vcf_G[[i]], MANE_G)
  query_hits <- queryHits(res)
  subject_hits <- subjectHits(res)
#Armar el data frame con los exones anotados, manteniendo la posicion dada por vcf  
result <- list()
  result[[i]] <- data.frame(
    chr = vcf[[i]]$CHROM[query_hits],
    start = vcf[[i]]$POS[query_hits],
    end = vcf[[i]]$END[query_hits],
    gen = MANE_file$gene_symbol[subject_hits],
    exon = MANE_file$exon[subject_hits]
  )
#Unir al data frame original los exones para no perder los calls en las regiones no anotadas por MANE y filtrar por lista de genes
final <- vcf[[i]] |>
  left_join(result[[i]], by = c("CHROM" = "chr", "POS" = "start", "END" = "end")) |>
  dplyr::select(CHROM, POS, END, ID, REF, ALT, QUAL, FILTER, GT, CN, NP, QA, QS, QSE, QSS, gen, exon) |>
  dplyr::filter(gen %in% genes_list)
#Guardar el resultado final en un archivo de texto
output_file[i] <- sub(".vcf", ".txt", filenames[i])
  print("Generado el anotado del VCF:")
  write.table(as.data.frame(final), file = output_file[i] ,quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)
}

print("END exon_annotation.R") 