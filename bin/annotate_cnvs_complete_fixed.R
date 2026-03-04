#!/usr/bin/env Rscript

# ============================================================
# ANÁLISIS COMPLETO DE CNVs SOMÁTICOS CON ANOTACIÓN DE GENES
# ============================================================
# Autor: [Tu nombre]
# Fecha: 2025-10-13
# Descripción: Filtra CNVs/LOH y anota genes de cáncer relevantes

library(optparse)
library(dplyr)
library(tidyr)  
library(GenomicRanges)

# ============================================================
# PARÁMETROS
# ============================================================
option_list <- list(
  make_option(c("--input_dir"), type="character", default="model_segments_results",
              help="Directorio con archivos .modelFinal.seg"),
  make_option(c("--mane"), type="character", default="~/anotado_somatic_CNV/exons_mane.txt",
              help="Archivo MANE_Select"),
  make_option(c("--output_dir"), type="character", default="CNV_annotated_results",
              help="Directorio de salida"),
  make_option(c("--min_size_kb"), type="numeric", default=50,
              help="Tamaño mínimo de segmento (kb)"),
  make_option(c("--maf_threshold"), type="numeric", default=0.35,
              help="Umbral para LOH (MinorAF < threshold)"),
  make_option(c("--log2cr_amp"), type="numeric", default=0.3,
              help="Umbral para amplificación (Log2CR > threshold)"),
  make_option(c("--log2cr_del"), type="numeric", default=-0.3,
              help="Umbral para deleción (Log2CR < threshold)"),
  make_option(c("--diploid_cn"), type="numeric", default=2,
              help="Número de copias diploide normal (default=2)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Crear directorio de salida
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# CARGAR DATOS
# ============================================================
cat("================================================\n")
cat("ANÁLISIS DE CNVs SOMÁTICOS CON ANOTACIÓN\n")
cat("================================================\n\n")

# Buscar archivos .modelFinal.seg
input_files <- list.files(opt$input_dir, 
                          pattern = ".*\\.modelFinal\\.seg$", 
                          full.names = TRUE)

if (length(input_files) == 0) {
  stop("No se encontraron archivos .modelFinal.seg en: ", opt$input_dir)
}

cat("Archivos encontrados:\n")
for (f in input_files) {
  cat("  -", basename(f), "\n")
}
cat("\n")

# Cargar MANE_select
cat("Cargando anotaciones MANE...\n")
MANE_file <- read.delim(opt$mane, sep = "\t", header = TRUE)
MANE_G <- makeGRangesFromDataFrame(
  MANE_file, 
  seqnames.field = "chr", 
  start.field = "start", 
  end.field = "end",
  keep.extra.columns = TRUE
)

# ============================================================
# GENES DE INTERÉS
# ============================================================
genes_list <- c(
  # Supresores tumorales
  "AIP", "ALK", "APC", "ATM", "ATR", "AXIN2", "BAP1", "BARD1", "BLM", 
  "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "BUB1B", "CASR", "CDC73", "CDH1", 
  "CDK4", "CDKN1B", "CDKN2A", "CFTR", "CHEK2", "CTNNA1", "CTRC", "DICER1", 
  "ENG", "EPCAM", "FANCC", "FANCD2", "FANCM", "FH", "FLCN", "GREM1", 
  "HOXB13", "KIF1B", "KIT", "LZTR1", "MAX", "MBD4", "MEN1", "MET", "MITF", 
  "MLH1", "MLH3", "MRE11A", "MSH2", "MSH3", "MSH6", "MUTYH", "NF1", "NF2", 
  "NBN", "NTHL1", "PALB2", "PDGFRA", "PIK3CA", "PMS2", "POLD1", "POLE", 
  "POT1", "PRKAR1A", "PRSS1", "PTCH1", "PTEN", "RAD51C", "RAD51D", "RB1", 
  "RECQL", "RECQL4", "RET", "RNF43", "RPS20", "SDHA", "SDHAF2", "SDHB", 
  "SDHC", "SDHD", "SMAD4", "SMARCA4", "SMARCB1", "SMARCE1", "SPINK1", 
  "STK11", "SUFU", "TMEM127", "TP53", "TSC1", "TSC2", "VHL", "WRN",
  # Oncogenes
  "ERBB2", "EGFR", "FGFR1", "FGFR2", "FGFR3", "KRAS", "BRAF", "MYC",
  "CCND1", "CCND2", "CCND3", "CCNE1", "CDK6", "MDM2", "MDM4", "AR"
)

# ============================================================
# FUNCIÓN PARA CLASIFICAR EVENTOS BASADA EN CN ABSOLUTO
# ============================================================
# Clasifica eventos usando CN_ABSOLUTE (corregido por pureza)
# en lugar de log2CR observado
classify_event_by_cn <- function(cn_absolute, minor_af, 
                                  diploid_cn = 2, maf_th = 0.35) {
  event <- case_when(
    # Deleción con LOH
    cn_absolute < diploid_cn & minor_af < maf_th ~ "DEL+LOH",
    
    # Deleción sin LOH (balanceada o parcial)
    cn_absolute < diploid_cn & minor_af >= maf_th ~ "DELETION",
    
    # Amplificación con LOH
    cn_absolute > diploid_cn & minor_af < maf_th ~ "AMP+LOH",
    
    # Amplificación sin LOH
    cn_absolute > diploid_cn & minor_af >= maf_th ~ "AMPLIFICATION",
    
    # Copy-neutral LOH (CN = 2 pero con LOH)
    cn_absolute == diploid_cn & minor_af < maf_th ~ "CN-LOH",
    
    # Normal (CN = 2 sin LOH)
    cn_absolute == diploid_cn & minor_af >= maf_th ~ "NEUTRAL",
    
    # Casos inesperados (no debería ocurrir)
    TRUE ~ "UNKNOWN"
  )
  return(event)
}

# ============================================================
# FUNCIÓN PARA ESTIMAR PUREZA TUMORAL
# ============================================================
# Método basado en Minor Allele Fraction (MAF) de segmentos con LOH
# Solo segmentos con LOH son informativos para pureza
# Pureza = 1 - (2 × MAF_promedio)
estimate_purity <- function(seg_data, min_segments = 1) {
  
  # Primero clasificar todos los segmentos para identificar LOH
  seg_classified <- seg_data %>%
    mutate(
      LOG2CR = LOG2_COPY_RATIO_POSTERIOR_50,
      MINOR_AF = MINOR_ALLELE_FRACTION_POSTERIOR_50,
      SIZE_KB = (END - START) / 1000
    )
  
  # Clasificar eventos temporalmente para filtrar LOH
  seg_classified <- seg_classified %>%
    mutate(
      HAS_LOH = case_when(
        is.na(MINOR_AF) ~ FALSE,
        MINOR_AF < 0.35 ~ TRUE,  # Umbral para LOH
        TRUE ~ FALSE
      )
    )
  
  # Filtrar solo segmentos con LOH
  loh_segments <- seg_classified %>%
    filter(
      HAS_LOH == TRUE,
      !is.na(MINOR_AF),
      MINOR_AF > 0,  # Valores válidos
      SIZE_KB >= 50  # Segmentos grandes para confiabilidad
    )
  
  if (nrow(loh_segments) < min_segments) {
    warning("No LOH segments detected for purity estimation (n=", nrow(loh_segments), ")")
    return(list(
      purity = NA, 
      n_segments = 0, 
      avg_maf = NA,
      weighted_avg_maf = NA,
      method = "no_LOH_detected"
    ))
  }
  
  # Calcular MAF promedio ponderado por tamaño del segmento
  loh_segments <- loh_segments %>%
    mutate(WEIGHT = SIZE_KB)
  
  weighted_avg_maf <- weighted.mean(
    loh_segments$MINOR_AF,
    loh_segments$WEIGHT
  )
  
  # Calcular MAF promedio simple
  avg_maf <- mean(loh_segments$MINOR_AF)
  
  # Fórmula: Pureza = 1 - (2 × MAF_promedio)
  # En tejido normal heterocigoto: MAF = 0.5
  # En tumor puro con LOH: MAF → 0
  # En mezcla: MAF = (1-p)*0.5 + p*0 = 0.5*(1-p)
  # Despejando: p = 1 - 2*MAF
  
  purity_estimate <- 1 - (2 * weighted_avg_maf)
  
  # Asegurar valores razonables (entre 0.1 y 1.0)
  # No forzar mínimo de 0.3 si los datos indican baja pureza
  purity_estimate <- max(0.1, min(1.0, purity_estimate))
  
  return(list(
    purity = round(purity_estimate, 4),
    n_segments = nrow(loh_segments),
    avg_maf = round(avg_maf, 4),
    weighted_avg_maf = round(weighted_avg_maf, 4),
    method = "LOH_based"
  ))
}

# ============================================================
# FUNCIÓN PARA CALCULAR NÚMERO DE COPIAS ABSOLUTO
# ============================================================
calculate_absolute_cn <- function(log2cr, purity, diploid_cn = 2) {
  
  # Asegurar que purity se recicle correctamente para el largo de log2cr
  if (length(purity) == 1) {
    purity <- rep(purity, length(log2cr))
  }
  
  # CN observado (mezcla tumor + normal)
  # Ratio observado = 2^(log2CR)
  ratio_observed <- 2^log2cr
  cn_observed <- diploid_cn * ratio_observed
  
  # Si NO hay pureza (sin LOH), devolver CN observado SIN corregir
  # Solo corregir cuando hay pureza válida
  cn_tumor <- ifelse(
    is.na(purity) | purity < 0.2 | purity > 1.0,
    cn_observed,  # Sin corrección: CN = 2 × 2^log2CR
    # Con corrección por pureza:
    (cn_observed - (1 - purity) * diploid_cn) / purity
  )
  
  # Redondear al entero más cercano
  cn_absolute <- round(cn_tumor)
  
  # Asegurar valores razonables (0-20)
  cn_absolute <- pmax(0, pmin(20, cn_absolute))
  
  return(cn_absolute)
}

# ============================================================
# PROCESAR CADA MUESTRA
# ============================================================
all_results <- list()
sample_names <- c()
purity_results <- list()

for (i in seq_along(input_files)) {
  file_path <- input_files[i]
  sample_name <- sub(".modelFinal.seg", "", basename(file_path))
  sample_names <- c(sample_names, sample_name)
  
  cat("Procesando:", sample_name, "\n")
  
  # Leer archivo
  seg_data <- read.delim(
    file_path,  
    sep = "\t", 
    header = TRUE, 
    comment.char = "@",
    stringsAsFactors = FALSE
  )
  
  # Estimar pureza tumoral
  cat("  - Estimando pureza tumoral...\n")
  purity_info <- estimate_purity(seg_data)
  sample_purity <- purity_info$purity
  
  # Guardar información de pureza
  purity_results[[i]] <- data.frame(
    SAMPLE = sample_name,
    PURITY = ifelse(is.na(sample_purity), NA, sample_purity),
    PURITY_PERCENT = ifelse(is.na(sample_purity), NA, round(sample_purity * 100, 2)),
    N_LOH_SEGMENTS = purity_info$n_segments,
    AVG_MAF = purity_info$avg_maf,
    WEIGHTED_AVG_MAF = purity_info$weighted_avg_maf,
    METHOD = purity_info$method
  )
  
  if (is.na(sample_purity)) {
    cat("    ⚠️  No se pudo calcular pureza (sin regiones LOH detectables)\n")
    cat("    CN se calculará SIN corrección por pureza (CN observado)\n")
    # NO asignar default - dejar como NA para que calculate_absolute_cn
    # devuelva CN observado sin corrección
  } else {
    cat("    Pureza estimada:", round(sample_purity * 100, 2), "%\n")
    cat("    Segmentos LOH analizados:", purity_info$n_segments, "\n")
  }
  
  # Filtrar y clasificar
  # IMPORTANTE: Primero calculamos CN_ABSOLUTE, luego clasificamos eventos
  # basándonos en el CN corregido por pureza (no en log2CR observado)
  seg_filtered <- seg_data %>%
    mutate(
      SIZE_KB = (END - START) / 1000,
      LOG2CR = LOG2_COPY_RATIO_POSTERIOR_50,
      MINOR_AF = MINOR_ALLELE_FRACTION_POSTERIOR_50,
      PURITY = sample_purity,
      # Primero: Calcular CN absoluto corregido por pureza
      CN_ABSOLUTE = calculate_absolute_cn(LOG2CR, sample_purity, opt$diploid_cn),
      # Segundo: Clasificar evento basándose en CN_ABSOLUTE
      EVENT_TYPE = classify_event_by_cn(CN_ABSOLUTE, MINOR_AF, opt$diploid_cn, opt$maf_threshold),
      SAMPLE = sample_name
    ) %>%
    filter(
      SIZE_KB >= opt$min_size_kb,
      EVENT_TYPE != "NEUTRAL"
    )
  
  # Convertir a GRanges
  seg_G <- makeGRangesFromDataFrame(
    seg_filtered, 
    seqnames.field = "CONTIG", 
    start.field = "START", 
    end.field = "END",
    keep.extra.columns = TRUE
  )
  
  # Anotar con MANE
  overlaps <- findOverlaps(seg_G, MANE_G)
  
  if (length(overlaps) > 0) {
    query_hits <- queryHits(overlaps)
    subject_hits <- subjectHits(overlaps)
    
    # Crear data frame anotado
    annotated <- seg_filtered[query_hits, ] %>%
      mutate(
        GENE = MANE_file$gene_symbol[subject_hits],
        EXON = MANE_file$exon_number[subject_hits]
      ) %>%
      filter(GENE %in% genes_list) %>%
      distinct() %>%
      arrange(CONTIG, START)
    
    all_results[[i]] <- annotated
  } else {
    cat("  ⚠ No se encontraron overlaps con genes\n")
    all_results[[i]] <- data.frame()
  }
}

# ============================================================
# COMBINAR RESULTADOS
# ============================================================
if (length(all_results) > 0) {
  combined_results <- bind_rows(all_results)
  
  # Guardar resultado completo
  output_combined <- file.path(opt$output_dir, "CNVs_annotated_all_samples.txt")
  write.table(
    combined_results, 
    file = output_combined,
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE, 
    col.names = TRUE
  )
  cat("\n✓ Resultados combinados guardados en:", output_combined, "\n")
  
  # ============================================================
  # GUARDAR RESULTADOS DE PUREZA TUMORAL
  # ============================================================
  purity_table <- bind_rows(purity_results)
  output_purity <- file.path(opt$output_dir, "tumor_purity_estimates.txt")
  write.table(
    purity_table,
    file = output_purity,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("✓ Estimaciones de pureza tumoral guardadas en:", output_purity, "\n")
  
  # ============================================================
  # GENERAR RESUMEN POR MUESTRA
  # ============================================================
  summary_by_sample <- combined_results %>%
    group_by(SAMPLE, EVENT_TYPE) %>%
    summarise(
      N_EVENTS = n(),
      N_GENES = n_distinct(GENE),
      TOTAL_SIZE_MB = sum(SIZE_KB) / 1000,
      .groups = "drop"
    ) %>%
    arrange(SAMPLE, EVENT_TYPE)
  
  output_summary <- file.path(opt$output_dir, "summary_by_sample.txt")
  write.table(
    summary_by_sample,
    file = output_summary,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("✓ Resumen por muestra guardado en:", output_summary, "\n")
  
  # ============================================================
  # GENERAR MATRIZ DE GENES x MUESTRAS
  # ============================================================
  gene_matrix <- combined_results %>%
    select(SAMPLE, GENE, EVENT_TYPE) %>%
    distinct() %>%
    pivot_wider(
      names_from = SAMPLE,
      values_from = EVENT_TYPE,
      values_fill = list(EVENT_TYPE = "NORMAL"),
      values_fn = function(x) paste(unique(x), collapse = ";")
    ) %>%
    arrange(GENE)
  
  output_matrix <- file.path(opt$output_dir, "gene_matrix.txt")
  write.table(
    gene_matrix,
    file = output_matrix,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("✓ Matriz gen x muestra guardada en:", output_matrix, "\n")
  
  # ============================================================
  # RESUMEN DE GENES MÁS AFECTADOS
  # ============================================================
  top_genes <- combined_results %>%
    group_by(GENE) %>%
    summarise(
      N_SAMPLES_AFFECTED = n_distinct(SAMPLE),
      EVENTS = paste(unique(EVENT_TYPE), collapse = ", "),
      AVG_LOG2CR = mean(LOG2CR),
      AVG_MINOR_AF = mean(MINOR_AF),
      .groups = "drop"
    ) %>%
    arrange(desc(N_SAMPLES_AFFECTED), GENE)
  
  output_top_genes <- file.path(opt$output_dir, "top_affected_genes.txt")
  write.table(
    top_genes,
    file = output_top_genes,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  cat("✓ Genes más afectados guardados en:", output_top_genes, "\n")
  
  # ============================================================
  # GUARDAR RESULTADOS INDIVIDUALES POR MUESTRA
  # ============================================================
  for (sample in unique(combined_results$SAMPLE)) {
    sample_data <- combined_results %>%
      filter(SAMPLE == sample) %>%
      select(-SAMPLE)
    
    output_individual <- file.path(opt$output_dir, paste0(sample, "_annotated.txt"))
    write.table(
      sample_data,
      file = output_individual,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  }
  cat("✓ Archivos individuales por muestra guardados\n")
  
  # ============================================================
  # IMPRIMIR RESUMEN EN CONSOLA
  # ============================================================
  cat("\n================================================\n")
  cat("RESUMEN DEL ANÁLISIS\n")
  cat("================================================\n\n")
  
  cat("Muestras procesadas:", length(sample_names), "\n")
  cat("Genes únicos afectados:", n_distinct(combined_results$GENE), "\n")
  cat("Total de eventos anotados:", nrow(combined_results), "\n\n")
  
  cat("Pureza tumoral estimada por muestra:\n")
  purity_display <- purity_table %>%
    select(SAMPLE, PURITY_PERCENT, N_LOH_SEGMENTS, WEIGHTED_AVG_MAF, METHOD)
  print(purity_display, row.names = FALSE)
  cat("\n")
  
  cat("Distribución de eventos:\n")
  print(table(combined_results$EVENT_TYPE))
  
  cat("\nDistribución de número de copias absoluto:\n")
  print(table(combined_results$CN_ABSOLUTE))
  
  cat("\nTop 10 genes más afectados:\n")
  print(head(top_genes, 10))
  
} else {
  cat("\n⚠ No se encontraron eventos que anotar\n")
}

cat("\n================================================\n")
cat("ANÁLISIS COMPLETADO\n")
cat("================================================\n")