# 필수 패키지
library(dplyr)
library(readr)

# 1. 파일 불러오기
deg1_ma <- read.csv("DEG1_Disease_High_vs_Low.csv")
deg1_sc <- read.csv("scRNA_DEG1_Disease_High_vs_Low.csv")

deg2_ma <- read.csv("DEG2_NonDisease_High_vs_Low.csv")
deg2_sc <- read.csv("scRNA_DEG2_NonDisease_High_vs_Low.csv")

deg3_ma <- read.csv("DEG3_High_Disease_vs_NonDisease.csv")
deg3_sc <- read.csv("scRNA_DEG3_High_Disease_vs_NonDisease.csv")

deg4_ma <- read.csv("DEG4_Low_Disease_vs_NonDisease.csv")
deg4_sc <- read.csv("scRNA_DEG4_Low_Disease_vs_NonDisease.csv")

# 2. 기준 적용 및 gene symbol 추출
get_filtered_genes <- function(df, fc_thresh = 1, p_thresh = 0.05) {
  # 유전자 심볼 열 추정
  gene_col <- intersect(c("hgnc_symbol", "gene", "symbol", "genes"), colnames(df))[1]
  # FC 및 p-value 열 추정
  fc_col   <- intersect(c("logFC", "avg_log2FC"), colnames(df))[1]
  p_col    <- intersect(c("adj.P.Val", "p_val_adj"), colnames(df))[1]
  
  if (is.null(gene_col) || is.null(fc_col) || is.null(p_col)) return(character(0))
  
  df %>%
    filter(abs(.data[[fc_col]]) > fc_thresh, .data[[p_col]] < p_thresh) %>%
    pull(.data[[gene_col]]) %>%
    tolower() %>%
    unique()
}



genes_deg1_ma <- get_filtered_genes(deg1_ma)
genes_deg1_sc <- get_filtered_genes(deg1_sc)
codeg1 <- intersect(genes_deg1_ma, genes_deg1_sc)

genes_deg2_ma <- get_filtered_genes(deg2_ma)
genes_deg2_sc <- get_filtered_genes(deg2_sc)
codeg2 <- intersect(genes_deg2_ma, genes_deg2_sc)

genes_deg3_ma <- get_filtered_genes(deg3_ma)
genes_deg3_sc <- get_filtered_genes(deg3_sc)
codeg3 <- intersect(genes_deg3_ma, genes_deg3_sc)

genes_deg4_ma <- get_filtered_genes(deg4_ma)
genes_deg4_sc <- get_filtered_genes(deg4_sc)
codeg4 <- intersect(genes_deg4_ma, genes_deg4_sc)


# 3. 교집합 (Co-DEG) 추출
codeg1 <- intersect(genes_deg1_ma, genes_deg1_sc)
codeg2 <- intersect(genes_deg2_ma, genes_deg2_sc)
codeg3 <- intersect(genes_deg3_ma, genes_deg3_sc)
codeg4 <- intersect(genes_deg4_ma, genes_deg4_sc)

# 4. 저장
write.csv(codeg1, "CoDEG1_Disease_High_vs_Low.csv", row.names = FALSE)
write.csv(codeg2, "CoDEG2_NonDisease_High_vs_Low.csv", row.names = FALSE)
write.csv(codeg3, "CoDEG3_High_Disease_vs_NonDisease.csv", row.names = FALSE)
write.csv(codeg4, "CoDEG4_Low_Disease_vs_NonDisease.csv", row.names = FALSE)



# 패키지
library(dplyr)
library(readr)

# 자동 병합 함수 정의
merge_codeg_info <- function(deg_ma, deg_sc,
                             fc_ma_col = "logFC", p_ma_col = "adj.P.Val", gene_ma_col = "hgnc_symbol",
                             fc_sc_col = "avg_log2FC", p_sc_col = "p_val_adj", gene_sc_col = "gene",
                             fc_thresh = 1, p_thresh = 0.05) {
  
  # 소문자 gene 컬럼 만들기
  deg_ma <- deg_ma %>%
    mutate(gene = tolower(.data[[gene_ma_col]])) %>%
    filter(abs(.data[[fc_ma_col]]) > fc_thresh, .data[[p_ma_col]] < p_thresh)
  
  deg_sc <- deg_sc %>%
    mutate(gene = tolower(.data[[gene_sc_col]])) %>%
    filter(abs(.data[[fc_sc_col]]) > fc_thresh, .data[[p_sc_col]] < p_thresh)
  
  # 공통 유전자 추출
  common_genes <- intersect(deg_ma$gene, deg_sc$gene)
  
  # 필요한 열만 추출
  deg_ma_sel <- deg_ma %>%
    filter(gene %in% common_genes) %>%
    select(gene, !!fc_ma_col, !!p_ma_col) %>%
    rename(logFC_microarray = !!fc_ma_col,
           adjP_microarray = !!p_ma_col)
  
  deg_sc_sel <- deg_sc %>%
    filter(gene %in% common_genes) %>%
    select(gene, !!fc_sc_col, !!p_sc_col) %>%
    rename(logFC_scRNA = !!fc_sc_col,
           adjP_scRNA = !!p_sc_col)
  
  # 병합 및 방향성 판단
  merged <- left_join(deg_ma_sel, deg_sc_sel, by = "gene") %>%
    mutate(direction = case_when(
      logFC_microarray > 0 & logFC_scRNA > 0 ~ "UpUp",
      logFC_microarray < 0 & logFC_scRNA < 0 ~ "DownDown",
      TRUE ~ "Discordant"
    ))
  
  return(merged)
}

# 각 파일 불러오기
deg1_ma <- read_csv("DEG1_Disease_High_vs_Low.csv")
deg1_sc <- read_csv("scRNA_DEG1_Disease_High_vs_Low.csv")

deg2_ma <- read_csv("DEG2_NonDisease_High_vs_Low.csv")
deg2_sc <- read_csv("scRNA_DEG2_NonDisease_High_vs_Low.csv")

deg3_ma <- read_csv("DEG3_High_Disease_vs_NonDisease.csv")
deg3_sc <- read_csv("scRNA_DEG3_High_Disease_vs_NonDisease.csv")

deg4_ma <- read_csv("DEG4_Low_Disease_vs_NonDisease.csv")
deg4_sc <- read_csv("scRNA_DEG4_Low_Disease_vs_NonDisease.csv")

# 병합 실행
codeg1_df <- merge_codeg_info(deg1_ma, deg1_sc)
codeg2_df <- merge_codeg_info(deg2_ma, deg2_sc)
codeg3_df <- merge_codeg_info(deg3_ma, deg3_sc)
codeg4_df <- merge_codeg_info(deg4_ma, deg4_sc)

# 저장
write.csv(codeg1_df, "CoDEG1_Annotated.csv", row.names = FALSE)
write.csv(codeg2_df, "CoDEG2_Annotated.csv", row.names = FALSE)
write.csv(codeg3_df, "CoDEG3_Annotated.csv", row.names = FALSE)
write.csv(codeg4_df, "CoDEG4_Annotated.csv", row.names = FALSE)







# ---------------------
# GO/KEGG 분석을 위해 Discordant gene 제외하고 UpUp + DownDown 유전자만 활용하도록 수정
# ---------------------
library(readr)
library(dplyr)
library(tools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(forcats)
library(patchwork)

# CoDEG 파일 목록
file_list <- list(
  "CoDEG1_Annotated.csv" = "CoDEG1: Disease High vs Low",
  "CoDEG2_Annotated.csv" = "CoDEG2: NonDisease High vs Low",
  "CoDEG3_Annotated.csv" = "CoDEG3: High Disease vs NonDisease",
  "CoDEG4_Annotated.csv" = "CoDEG4: Low Disease vs NonDisease"
)

# 유전자 리스트 불러오기 함수 (UpUp/DownDown만)
load_gene_symbols <- function(df) {
  gene_symbols <- df %>%
    filter(direction %in% c("UpUp", "DownDown")) %>%
    pull(gene) %>%
    na.omit() %>% trimws() %>% unique() %>% .[. != ""] %>% toupper()
  valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  gene_symbols[gene_symbols %in% valid_symbols]
}

# GO 분석 함수
run_go_three_plot <- function(gene_symbols, condition_name, title_text = "GO Enrichment") {
  # SYMBOL → ENTREZID 변환
  gene_entrez <- suppressMessages(
    bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
      distinct(ENTREZID) %>%
      pull(ENTREZID)
  )
  
  if (length(gene_entrez) < 5) {
    warning(paste(condition_name, "유효한 유전자 수 부족 → 건너뜀"))
    return(NULL)
  }
  
  # GO 서브함수 (6개로 제한)
  get_go_df <- function(ont) {
    go <- enrichGO(gene_entrez, OrgDb = org.Hs.eg.db, ont = ont,
                   pAdjustMethod = "BH", readable = TRUE)
    if (is.null(go) || nrow(go@result) == 0) return(NULL)
    go@result %>%
      slice_min(p.adjust, n = 6) %>%
      mutate(
        GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) /
          as.numeric(sub(".*/", "", GeneRatio)),
        Description = fct_reorder(Description, GeneRatio)
      )
  }
  
  df_bp <- get_go_df("BP")
  df_cc <- get_go_df("CC")
  df_mf <- get_go_df("MF")
  
  # Barplot 생성 함수
  make_barplot <- function(df, go_title) {
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x = Description, y = GeneRatio)) +
      geom_col(aes(fill = -log10(p.adjust)), width = 0.7, color = "black") +
      geom_line(aes(y = -log10(p.adjust) / max(-log10(p.adjust)) * max(GeneRatio), group = 1),
                color = "#FFD700", linewidth = 1.2) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.6) +
      coord_flip(clip = "off") +
      scale_fill_gradient(low = "#92c5de", high = "#ca0020") +
      labs(title = go_title, y = "GeneRatio + scaled -log10(p)", x = NULL) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid = element_line(color = "grey90"),
        axis.text.y = element_text(size = 12, face = "bold")  # 축 텍스트 강조
      )
  }
  
  # 최종 플롯 조합
  final_plot <- (make_barplot(df_bp, "Biological Process") /
                   make_barplot(df_cc, "Cellular Component") /
                   make_barplot(df_mf, "Molecular Function")) +
    plot_annotation(title = title_text) &
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  
  # 저장
  ggsave(paste0("GO_Enrichment_", condition_name, ".png"),
         plot = final_plot, width = 10, height = 12, dpi = 600)
  cat("✅ GO 저장 완료:", paste0("GO_Enrichment_", condition_name, ".png"), "\n")
}




# KEGG 분석 함수
  run_kegg_plot <- function(gene_symbols, condition_name, organism = "hsa") {
    # 1. SYMBOL → ENTREZID 변환
    gene_entrez <- suppressMessages(
      bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
        distinct(ENTREZID) %>%
        pull(ENTREZID)
    )
    
    if (length(gene_entrez) < 5) {
      warning(paste(condition_name, "유효한 유전자 수 부족 → 건너뜀"))
      return(NULL)
    }
    
    # 2. KEGG 분석
    kegg <- enrichKEGG(gene = gene_entrez, organism = organism, pAdjustMethod = "BH")
    
    if (is.null(kegg) || nrow(kegg@result) == 0) {
      warning(paste(condition_name, "KEGG 결과 없음"))
      return(NULL)
    }
    
    # 3. 상위 15개 추출 및 전처리
    df <- kegg@result %>%
      slice_min(p.adjust, n = 15) %>%
      mutate(
        GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) /
          as.numeric(sub(".*/", "", GeneRatio)),
        Description = fct_reorder(Description, GeneRatio)
      )
    
    # 4. 플롯 생성
    p <- ggplot(df, aes(x = Description, y = GeneRatio)) +
      geom_col(aes(fill = -log10(p.adjust)), width = 0.7, color = "black") +
      geom_line(aes(y = -log10(p.adjust) / max(-log10(p.adjust)) * max(GeneRatio), group = 1),
                color = "#FFD700", linewidth = 1.2) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.6) +
      coord_flip(clip = "off") +
      scale_fill_gradient(low = "#92c5de", high = "#ca0020") +
      labs(
        title = paste("KEGG Pathway Enrichment:", condition_name),
        y = "GeneRatio + scaled -log10(p)", x = NULL
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid = element_line(color = "grey90"),
        axis.text.y = element_text(size = 12, face = "bold")
      )
    
    # 5. 저장
    ggsave(paste0("KEGG_Enrichment_", condition_name, ".png"),
           plot = p, width = 10, height = 8, dpi = 600)
    cat("KEGG 저장 완료:", paste0("KEGG_Enrichment_", condition_name, ".png"), "\n")
  }
  

# 실행: UpUp + DownDown만 사용하여 GO/KEGG
for (f in names(file_list)) {
  df <- read_csv(f, show_col_types = FALSE)
  gene_symbols <- load_gene_symbols(df)
  condition <- file_path_sans_ext(basename(f))
  run_go_three_plot(gene_symbols, condition, paste("GO Enrichment -", file_list[[f]]))
  run_kegg_plot(gene_symbols, condition)
}

# ─────────────────────────────────────────────
# 1) 패키지 로드
# ─────────────────────────────────────────────
library(tidyverse)
library(fmsb)
library(scales)
library(stringr)

# 2) 데이터 불러오기  (−log10 p) 변환
kegg_df <- read_csv("kegg_pvals.csv", show_col_types = FALSE) %>%
  mutate(score = -log10(pvalue))

# 3) 조건별 top‑10 경로 wide 변환
top_wide <- kegg_df %>%
  group_by(CoDEG) %>%
  slice_min(pvalue, n = 10, with_ties = FALSE) %>%  # ← 10개 유지
  ungroup() %>%
  select(CoDEG, Pathway, score) %>%
  mutate(Pathway = str_replace_all(Pathway, "[^A-Za-z0-9]+", "_")) %>%
  pivot_wider(names_from = Pathway, values_from = score, values_fill = 0) %>%
  arrange(CoDEG)

# 4) 0‑1 스케일링 + min/max 행
score_cols <- setdiff(names(top_wide), "CoDEG")
top_scaled <- top_wide
top_scaled[score_cols] <- lapply(top_scaled[score_cols],
                                 \(x) (x - min(x)) / (max(x) - min(x) + 1e-12))

radar_mat   <- top_scaled %>% select(all_of(score_cols)) %>% as.data.frame()
radar_ready <- rbind(rep(1, ncol(radar_mat)),
                     rep(0, ncol(radar_mat)),
                     radar_mat)
rownames(radar_ready) <- c("max", "min", paste0("CoDEG", top_scaled$CoDEG))

# 5) 색상 + 라벨(줄바꿈 24자) + 폰트 크기
colors_border <- c("firebrick", "darkorange", "steelblue", "cyan4")
colors_fill   <- alpha(colors_border, 0.25)
vlabels       <- gsub("_", " ", score_cols) |> str_wrap(24)

vl_size  <- 0.45   # 축 라벨
cal_size <- 0.45   # 100/75/50/25
leg_size <- 0.80   # 범례

# 6) 레이더플롯 & 저장
png("Radar_KEGG_CoDEG1-4.png", width = 2400, height = 2000, res = 300)
par(font = 2, mar = c(1,2,2,1))  # 볼드 + 여백

radarchart(radar_ready,
           axistype = 1,
           vlabels  = vlabels,
           vlcex    = vl_size,
           calcex   = cal_size,
           pcol     = colors_border,
           pfcol    = colors_fill,
           plwd     = 2,
           plty     = 1,
           cglcol   = "grey65",
           cglty    = 1,
           cglwd    = 0.7,
           title    = "KEGG spectrum (CoDEG1–4)\n(top 10 pathways)")

legend("topright", legend = paste0("CoDEG", 1:4),
       bty = "n", pch = 20, col = colors_border,
       text.col = "black", text.font = 1,
       cex = leg_size, pt.cex = 1.2)

dev.off()
cat("Radar_KEGG_CoDEG1-4.png 저장 완료 (top‑10 + 더 작은 글씨)\n")


# 패키지
library(readr)
library(dplyr)
library(fmsb)
library(stringr)
library(scales)   # alpha()

#  레이더 플롯 함수
draw_go_radar <- function(wide_csv,
                          title_prefix   = "GO–",
                          out_png        = NULL,
                          zigzag_labels  = TRUE,
                          vl_size        = 0.45,
                          cal_size       = 0.40,
                          leg_size       = 0.80) {
  
  ## 1. 파일 읽기 (CoDEG + 0–1 스케일 컬럼)
  wide_df <- read_csv(wide_csv, show_col_types = FALSE)
  
  ## 2. fmsb용 min/max 행 추가
  score_cols <- setdiff(names(wide_df), "CoDEG")
  radar_mat  <- wide_df %>% select(all_of(score_cols)) %>% as.data.frame()
  radar_ready <- rbind(rep(1, ncol(radar_mat)),  # max
                       rep(0, ncol(radar_mat)),  # min
                       radar_mat)
  rownames(radar_ready) <- c("max", "min", paste0("CoDEG", wide_df$CoDEG))
  
  ## 3. 레이블 처리 (언더바→공백, 줄바꿈)
  vlabels <- gsub("_", " ", score_cols) |> stringr::str_wrap(24)
  if (zigzag_labels) {
    vlabels <- ifelse(seq_along(vlabels) %% 2 == 0,
                      paste0("\n", vlabels), vlabels)
  }
  
  ## 4. 색상
  colors_border <- c("firebrick", "darkorange", "steelblue", "cyan4")[seq_along(wide_df$CoDEG)]
  colors_fill   <- alpha(colors_border, 0.25)
  
  ## 5. 출력 파일 이름
  if (is.null(out_png)) {
    onto <- stringr::str_extract(wide_csv, "(BP|CC|MF)")
    out_png <- paste0("Radar_", title_prefix, onto, ".png")
  }
  
  ## 6. 그리기 & 저장
  png(out_png, width = 2400, height = 2000, res = 300)
  par(font = 2, mar = c(1, 2, 2, 1))
  
  radarchart(radar_ready,
             axistype = 1,
             vlabels  = vlabels,
             vlcex    = vl_size,
             calcex   = cal_size,
             pcol     = colors_border,
             pfcol    = colors_fill,
             plwd     = 2,
             plty     = 1,
             cglcol   = "grey65",
             cglty    = 1,
             cglwd    = 0.7,
             title    = paste0(title_prefix,
                               stringr::str_extract(wide_csv, "(BP|CC|MF)"),
                               " spectrum (CoDEG1–4)\n(top 10 pathways)"))
  
  legend("topright", legend = paste0("CoDEG", wide_df$CoDEG),
         bty = "n", pch = 20, col = colors_border,
         text.col = "black", text.font = 1,
         cex = leg_size, pt.cex = 1.2)
  
  dev.off()
  cat("저장 완료:", out_png, "\n")
}

#  호출 예시
draw_go_radar("GO_BP_wide.csv")   # 생물학적 과정
draw_go_radar("GO_CC_wide.csv")   # 세포 소구조
draw_go_radar("GO_MF_wide.csv")   # 분자 기능





# ---------------------
# Volcano Plot
# ---------------------
library(ggplot2)
library(readr)
library(dplyr)
library(tools)

plot_volcano <- function(file_path, title_label = NULL) {
  df <- read_csv(file_path, show_col_types = FALSE)
  
  # 컬럼명 표준화
  df <- df %>%
    rename(
      gene = gene,
      log2fc = logFC_microarray,
      pvalue = adjP_microarray,
      direction = direction
    )
  
  # 필수 컬럼 확인
  if (!all(c("gene", "log2fc", "pvalue", "direction") %in% colnames(df))) {
    warning(paste("필수 컬럼이 누락됨:", file_path))
    return(NULL)
  }
  
  # Volcano용 데이터 가공
  df <- df %>%
    mutate(
      log10_p = -log10(pvalue),
      sig = case_when(
        direction == "UpUp" ~ "Upregulated",
        direction == "DownDown" ~ "Downregulated",
        TRUE ~ "Discordant/Other"
      )
    )
  
  title_text <- ifelse(is.null(title_label), tools::file_path_sans_ext(basename(file_path)), title_label)
  
  # Volcano plot 생성
  p <- ggplot(df, aes(x = log2fc, y = log10_p, color = sig)) +
    # 배경 색상 geom_rect 추가
    annotate("rect", xmin = -Inf, xmax = -1, ymin = 0, ymax = Inf,
             alpha = 0.1, fill = "#0571b0") +
    annotate("rect", xmin = 1, xmax = Inf, ymin = 0, ymax = Inf,
             alpha = 0.1, fill = "#ca0020") +
    
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Upregulated" = "#ca0020", "Downregulated" = "#0571b0", "Discordant/Other" = "gray")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.6) +
    labs(
      title = paste("Volcano Plot -", title_text),
      x = "log2 Fold Change (Microarray)",
      y = "-log10(adj.P)",
      color = "Significance"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  # 파일 저장
  ggsave(paste0("Volcano_", tools::file_path_sans_ext(basename(file_path)), ".png"),
         plot = p, width = 8, height = 6, dpi = 400)
  
  cat("Volcano plot 저장 완료:", file_path, "\n")
}

# 반복 실행
file_list <- list(
  "CoDEG1_Annotated.csv" = "CoDEG1: Disease High vs Low",
  "CoDEG2_Annotated.csv" = "CoDEG2: NonDisease High vs Low",
  "CoDEG3_Annotated.csv" = "CoDEG3: High Disease vs NonDisease",
  "CoDEG4_Annotated.csv" = "CoDEG4: Low Disease vs NonDisease"
)

for (f in names(file_list)) {
  plot_volcano(f, title_label = file_list[[f]])
}





# ---------------------
# PPI 유사 네트워크 분석
# --------------------
# 필수 패키지 로드
library(dplyr)
library(purrr)
library(igraph)

# 1. 중심 유전자 목록 (염증 및 죽상동맥경화 관련)
key_genes <- c(
  "OASL", "CD14", "CD68", "CD3D", "CD3E", "CD8A", "CD27",
  "STAT1", "IFIT1", "IFI44L", "MX1", "IRF7", "IRF8",
  "HLA-DRA", "HLA-DPB1", "HLA-DQA1", "TYROBP", "LCP2",
  "TLR2", "NLRP3", "CCL5", "CXCR4", "CXCL10",
  "TNF", "PTPRC", "ITGAL", "ITGAM", "SPI1", "IL2RG", "CD86"
)

# 2. 각 유전자에서 상관계수 기준 상위 2개만 선택
subset_edges <- map_dfr(key_genes, function(g) {
  cor_df %>%
    filter(from == g | to == g) %>%
    arrange(desc(abs(correlation))) %>%
    slice_head(n = 2)
})

# 3. 네트워크 객체 생성
g <- graph_from_data_frame(subset_edges, directed = FALSE)

# 4. 허브 유전자 추출 (연결도 기준)
hub_scores <- degree(g, mode = "all")
hub_genes <- sort(hub_scores, decreasing = TRUE)
top_hub_genes <- head(hub_genes, 10)  # 상위 10개

# 5. 시각화 속성 설정
V(g)$color <- ifelse(V(g)$name %in% names(top_hub_genes), "tomato", "skyblue")
V(g)$size <- ifelse(V(g)$name %in% names(top_hub_genes), 25, 12)
V(g)$label.cex <- ifelse(V(g)$name %in% names(top_hub_genes), 1.6, 1.0)
V(g)$label.font <- ifelse(V(g)$name %in% names(top_hub_genes), 2, 1)  # 허브 유전자 굵게
V(g)$label.color <- "black"

E(g)$width <- abs(E(g)$correlation) * 4
E(g)$color <- "gray40"

# 6. 레이아웃 & 플롯
set.seed(123)
layout <- layout_with_fr(g)

plot(g,
     layout = layout,
     vertex.label.family = "sans",
     edge.curved = 0.1,
     main = "Gene Correlation Network with Highlighted Hub Genes")




# ---------------------
# GSVA
# --------------------
# 패키지 로드
library(GSVA)
library(msigdbr)
library(BiocParallel)

# 1. 발현 행렬 로딩
expr <- readr::read_csv("GSE100927_gene_mapped_expression.csv", show_col_types = FALSE) |>
  as.data.frame()
rownames(expr) <- toupper(expr$symbol)
expr <- expr[, -1]

# 2. Gene Set 준비 (MSigDB Hallmark)
msig <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(msig$gene_symbol, msig$gs_name)

# 3. GSVA 파라미터 객체 만들기 (※ geneSets 인자 사용)
gsvapar <- gsvaParam(expr = as.matrix(expr),
                     geneSets = gene_sets,
                     kcdf = "Gaussian")

# 4. GSVA 실행
gsva_result <- gsva(gsvapar)

# 5. 결과 확인
head(gsva_result[, 1:5])



# 1. OASL 발현 추출 및 그룹 지정 (수치형 변환 포함)
oasl_expr <- as.numeric(expr["OASL", ])
names(oasl_expr) <- colnames(expr)

# 그룹 지정
oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")
names(oasl_group) <- colnames(expr)



# 추가 패키지
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(tibble)

# 1. OASL 발현 추출 및 그룹 지정
oasl_expr <- expr["OASL", ]
oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")
names(oasl_group) <- colnames(expr)

# 2. Z-score 정규화 (pathway-wise)
gsva_z <- t(scale(t(gsva_result)))

# 3. 히트맵 (OASL High/Low annotation 포함)
column_ha <- HeatmapAnnotation(
  OASL = oasl_group[colnames(gsva_z)],
  col = list(OASL = c(High = "#E41A1C", Low = "#377EB8"))
)

Heatmap(gsva_z,
        name = "GSVA Z-score",
        top_annotation = column_ha,
        show_column_names = FALSE,
        show_row_names = TRUE,
        cluster_columns = TRUE,
        cluster_rows = TRUE,
        column_title = "GSVA Pathway Activity (Z-score)",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"))

# 4. Boxplot용 데이터 변환
gsva_df <- as.data.frame(t(gsva_result))
gsva_df$Sample <- rownames(gsva_df)
gsva_df$OASL_Group <- oasl_group[gsva_df$Sample]

gsva_long <- gsva_df %>%
  pivot_longer(-c(Sample, OASL_Group), names_to = "Pathway", values_to = "GSVA_Score")

# 5. 선택된 상위 항목 (예: 염증 관련 6개)
top_pathways <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT"
)

gsva_top <- gsva_long %>% filter(Pathway %in% top_pathways)

# 6. Boxplot 시각화
ggplot(gsva_top, aes(x = OASL_Group, y = GSVA_Score, fill = OASL_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~ Pathway, scales = "free_y") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("High" = "#E41A1C", "Low" = "#377EB8")) +
  labs(title = "Pathway Activity by OASL Expression Group (GSVA)",
       x = "OASL Expression Group", y = "GSVA Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


library(forcats)

# pathway × group 평균값 계산
gsva_summary <- gsva_df %>%
  group_by(Pathway, Group) %>%
  summarise(MeanScore = mean(Score), .groups = "drop")

# pathway 순서 지정: 평균 차이 큰 순
pathway_order <- gsva_summary %>%
  pivot_wider(names_from = Group, values_from = MeanScore) %>%
  mutate(diff = abs(High - Low)) %>%
  arrange(desc(diff)) %>%
  pull(Pathway)

# 막대그래프
ggplot(gsva_summary, aes(x = MeanScore, y = fct_inorder(Pathway), fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("High" = "firebrick", "Low" = "steelblue")) +
  labs(title = "GSVA Pathway Activity (Mean Score by OASL Group)",
       x = "Mean GSVA Score", y = "Pathway") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "top")


# 1. GSVA 결과 전치 (샘플 × 패스웨이)
gsva_mat <- t(gsva_result)  # [샘플 x 경로] 형태

# 2. Z-score 변환 (각 경로별로 z-score)
gsva_z <- scale(gsva_mat)  # scale은 열 단위로 처리하므로, 패스웨이별 Z-score로 적절함

# 3. OASL 발현 추출 및 그룹화
oasl_expr <- as.numeric(expr["OASL", ])
names(oasl_expr) <- colnames(expr)
oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")

# 4. 샘플 정렬 (High → Low)
ordered_samples <- names(sort(oasl_group, decreasing = TRUE))
gsva_z <- gsva_z[ordered_samples, ]  # 샘플 순서 적용

# 5. 주석 및 색상 정의
ann <- data.frame(Group = oasl_group[ordered_samples])
rownames(ann) <- ordered_samples
ann_colors <- list(Group = c(High = "firebrick", Low = "steelblue"))

# 6. 히트맵 시각화
pheatmap::pheatmap(t(gsva_z),  # [패스웨이 x 샘플]
                   annotation_col = ann,
                   annotation_colors = ann_colors,
                   cluster_cols = FALSE,  # 샘플 정렬 유지
                   cluster_rows = TRUE,
                   show_colnames = FALSE,
                   main = "GSVA Z-score Heatmap (OASL High → Low)",
                   fontsize_row = 8)


library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)
library(RColorBrewer)

# 염증 관련 pathway 지정
inflammatory_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_COAGULATION",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_TGF_BETA_SIGNALING"
)

# 데이터 필터링 및 평균 계산
gsva_summary <- gsva_df %>%
  filter(Pathway %in% inflammatory_pathways) %>%
  group_by(Pathway, Group) %>%
  summarise(MeanScore = mean(Score), .groups = "drop")

# 정렬: High vs Low 차이 기준
pathway_order <- gsva_summary %>%
  pivot_wider(names_from = Group, values_from = MeanScore) %>%
  mutate(diff = abs(High - Low)) %>%
  arrange(desc(diff)) %>%
  pull(Pathway)

# 색상 팔레트
group_colors <- c("High" = "#D73027", "Low" = "#4575B4")  # Red-blue 계열

# 시각화
ggplot(gsva_summary, aes(x = MeanScore, y = fct_relevel(Pathway, pathway_order), fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "GSVA Pathways (OASL High vs Low)",
    x = "Mean GSVA Score",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 15),
    axis.text.x = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


library(pheatmap)
library(dplyr)

# 염증 관련 pathway만 선택 + 정렬
gsva_sel <- gsva_result[inflammatory_pathways, ]

# rowname에서 "HALLMARK_" 제거
rownames(gsva_sel) <- gsub("^HALLMARK_", "", rownames(gsva_sel))

# Z-score scaling (pathway 별)
gsva_scaled <- t(scale(t(gsva_sel)))  # row 단위 scale 후 transpose

# OASL 기준 그룹 정의
oasl_expr <- expr["OASL", ]
oasl_expr <- as.numeric(oasl_expr)
oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")
names(oasl_group) <- colnames(expr)

# 샘플 정렬
ordered_samples <- names(sort(oasl_group, decreasing = TRUE))
gsva_scaled <- gsva_scaled[, ordered_samples]

# 어노테이션 및 색상
ann <- data.frame(Group = oasl_group[ordered_samples])
rownames(ann) <- ordered_samples
ann_colors <- list(Group = c(High = "firebrick", Low = "steelblue"))

# 히트맵 그리기
pheatmap(gsva_scaled,
         annotation_col = ann,
         annotation_colors = ann_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize_row = 12,
         color = colorRampPalette(c("#313695", "white", "#A50026"))(100),
         main = "Z-score Heatmap of Inflammatory Pathways (GSVA, OASL High → Low)",
         border_color = "grey70")



# ---------------------
# TF 분석
# --------------------
# 패키지 로드
library(readr)
library(dplyr)
library(viper)
library(dorothea)
library(pheatmap)
library(ggplot2)
library(forcats)
library(tidyr)

# 1. 발현 행렬 불러오기 --------------------------------------------------
expr <- read_csv("GSE100927_gene_mapped_expression.csv", show_col_types = FALSE) |> as.data.frame()
rownames(expr) <- toupper(expr$symbol)
expr <- expr[, -1]

# 2. regulon 로딩 및 변환 ------------------------------------------------
data(dorothea_hs, package = "dorothea")
regulon_df <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))

convert_regulon_to_viper <- function(dorothea_regulon) {
  regulon_split <- split(dorothea_regulon, dorothea_regulon$tf)
  regulon_list <- lapply(regulon_split, function(reg) {
    tfmode <- setNames(reg$mor, reg$target)
    likelihood <- setNames(abs(reg$mor), reg$target)
    list(tfmode = tfmode, likelihood = likelihood)
  })
  return(regulon_list)
}

regulon_list <- convert_regulon_to_viper(regulon_df)

# 3. TF activity 계산 -----------------------------------------------------
tf_activity <- viper(as.matrix(expr), regulon_list, method = "scale")

# 4. OASL 기반 그룹 분리 --------------------------------------------------
oasl_expr <- as.numeric(expr["OASL", ])
oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")
names(oasl_group) <- colnames(expr)

# TF activity 행렬 정렬
tf_z <- t(scale(t(tf_activity)))  # Z-score

# 샘플 순서 고정
ordered_samples <- names(sort(oasl_group, decreasing = TRUE))
tf_z <- tf_z[, ordered_samples]
annotation_col <- data.frame(Group = oasl_group[ordered_samples])
rownames(annotation_col) <- ordered_samples
ann_colors <- list(Group = c(High = "firebrick", Low = "steelblue"))

# 5. 시각화: 히트맵 ------------------------------------------------------
# 1. OASL 발현 추출 및 그룹 정의
oasl_expr <- as.numeric(expr["OASL", ])
names(oasl_expr) <- colnames(expr)
oasl_group <- ifelse(oasl_expr >= median(oasl_expr, na.rm = TRUE), "High", "Low")

# 2. 샘플 순서 지정: High → Low 순
ordered_samples <- names(sort(oasl_group, decreasing = TRUE))

# 3. TF activity Z-score 행렬 정렬
tf_z <- t(scale(t(tf_activity)))   # TF × Sample
tf_z <- tf_z[, ordered_samples]   # 샘플 정렬 적용

# 4. annotation 설정
annotation_col <- data.frame(Group = oasl_group[ordered_samples])
rownames(annotation_col) <- ordered_samples
ann_colors <- list(Group = c(High = "firebrick", Low = "steelblue"))

# 5. 히트맵 시각화
pheatmap(tf_z,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         show_rownames = TRUE,
         clustering_method = "ward.D2",  # TF만 클러스터링
         cluster_cols = FALSE,           # 샘플은 클러스터링 하지 않음
         main = "TF Activity Heatmap (Z-score, Ordered by OASL)",
         fontsize_row = 8,
         border_color = "gray90")

# 6. 시각화: 바플롯 (상위 변화 TF만 선택) -------------------------------
# DF로 변환
tf_df <- as.data.frame(t(tf_activity)) |> rownames_to_column("Sample")
tf_df$Group <- oasl_group[tf_df$Sample]
tf_df <- pivot_longer(tf_df, -c(Sample, Group), names_to = "TF", values_to = "Activity")

# 평균 차이 큰 TF top 15
tf_summary <- tf_df %>%
  group_by(TF, Group) %>%
  summarise(MeanScore = mean(Activity), .groups = "drop")

top_tfs <- tf_summary %>%
  pivot_wider(names_from = Group, values_from = MeanScore) %>%
  mutate(diff = abs(High - Low)) %>%
  arrange(desc(diff)) %>%
  slice_head(n = 15) %>%
  pull(TF)

tf_summary <- tf_summary %>% filter(TF %in% top_tfs)

# 바플롯
ggplot(tf_summary, aes(x = MeanScore, y = fct_reorder(TF, MeanScore), fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("High" = "firebrick", "Low" = "steelblue")) +
  labs(title = "Top 15 TF Activities by OASL Group", x = "Mean Activity", y = "Transcription Factor") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.text.y = element_text(face = "bold"),
        legend.position = "top",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))







# ---------------------
# 머신 러닝
# --------------------
# 필수 패키지 로드
library(readr)
library(dplyr)
library(glmnet)
library(caret)
library(randomForest)
library(pROC)

# 1. 발현 데이터 불러오기 --------------------------------------------------
expr <- read_csv("GSE100927_gene_mapped_expression.csv", show_col_types = FALSE)
expr <- as.data.frame(expr)
rownames(expr) <- toupper(expr$symbol)
expr <- expr[, -1]  # symbol 컬럼 제거

# 2. 그룹 지정 (Control 포함 여부로 Disease vs NonDisease 분리) ------------
group <- ifelse(grepl("Control", colnames(expr)), "NonDisease", "Disease")
group <- factor(group)
table(group)

# 3. CoDEG3 + OASL 유전자 리스트 불러오기 -----------------------------------
gene_list <- scan("CoDEG4_DownDown_genes.txt", what = character())
gene_list <- unique(toupper(gene_list))  # 대문자로 정규화

# 4. 발현 행렬 필터링 + 그룹 병합 ------------------------------------------
expr_sub <- expr[rownames(expr) %in% gene_list, ]
expr_sub <- t(expr_sub)
expr_sub <- as.data.frame(expr_sub)
expr_sub$group <- group

# 5. stratified train/test 분할 ---------------------------------------------
set.seed(123)
train_idx <- createDataPartition(expr_sub$group, p = 0.7, list = FALSE)
train <- expr_sub[train_idx, ]
test  <- expr_sub[-train_idx, ]

# 6. LASSO 모델 학습을 위한 x/y 준비 ----------------------------------------
x_train <- train[, !colnames(train) %in% "group"]
x_train <- x_train %>% mutate(across(everything(), as.numeric))
x_train <- x_train[, apply(x_train, 2, var, na.rm = TRUE) != 0]
x_train_matrix <- as.matrix(x_train)

y_train <- ifelse(train$group == "Disease", 1, 0)
x_df <- as.data.frame(x_train_matrix)
x_df$group <- factor(ifelse(y_train == 1, "Disease", "NonDisease"), levels = c("NonDisease", "Disease"))

# 7. caret 기반 LASSO 모델 학습 ---------------------------------------------
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

lasso_model <- train(
  group ~ ., data = x_df,
  method = "glmnet",
  family = "binomial",
  metric = "ROC",
  trControl = ctrl,
  tuneGrid = expand.grid(alpha = 1, lambda = seq(0.0001, 1, length = 100))
)

# 8. 최적 lambda와 유전자 선택 확인 -----------------------------------------
best_lambda <- lasso_model$bestTune$lambda
cat("선택된 lambda:", best_lambda, "\n")

coefs <- coef(lasso_model$finalModel, s = best_lambda)
lasso_genes <- rownames(coefs)[which(coefs != 0)]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")

cat("선택된 유전자 목록:\n")
print(lasso_genes)

# 9. Random Forest
# group은 factor여야 함
train$group <- factor(train$group)

# 특수문자 포함 변수명을 안전하게 사용하기 위해 수동 지정
x_train <- train[, !(colnames(train) == "group")]
y_train <- train$group

rf_model <- randomForest(x = x_train, y = y_train, importance = TRUE)

# 중요도 상위 유전자 추출
rf_importance <- importance(rf_model)
rf_genes <- rownames(rf_importance)[order(rf_importance[, 1], decreasing = TRUE)][1:20]

cat("RF 상위 유전자:\n")
print(rf_genes)

# 10. SVM-RFE
ctrl <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
svm_profile <- rfe(train[, -ncol(train)], train$group,
                   sizes = c(5, 10, 20), rfeControl = ctrl)
svm_genes <- predictors(svm_profile)
cat("SVM-RFE 유전자:\n")
print(svm_genes)

# 11. ROC 분석 (LASSO)
lasso_test_expr <- test[, lasso_genes, drop = FALSE]



# Heatmap
# 의미: 유전자 간의 상관 관계를 시각적으로 보여준다. 상관관계가 높은 유전자들은 색상이 비슷하며, 낮은 유전자들은 다른 색을 가진다.
# 활용: 여러 유전자 간 상관 관계를 시각화하여 유전자 그룹 또는 모듈을 식별할 수 있다.
library(pheatmap)

# 상위 20개 유전자 간의 상관 관계
corr_matrix <- cor(train[, lasso_genes])
pheatmap(corr_matrix, cluster_rows = TRUE, cluster_cols = TRUE, main = "유전자 간 상관 관계")


# Correlation Plot
# 의미: 유전자 간 상관 관계를 점으로 시각화하여 유전자 간의 연관성을 빠르게 확인할 수 있다.
# 활용: 특정 유전자들이 다른 유전자들과 얼마나 연관이 있는지 파악할 수 있다.
library(corrplot)

# 상위 20개 유전자 간의 상관 관계
corr_matrix <- cor(train[, lasso_genes])
corrplot(corr_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45)



# 두 모델에서 겹치는 유전자
common_lasso_rf <- intersect(lasso_genes, rf_genes)
common_lasso_svm <- intersect(lasso_genes, svm_genes)
common_rf_svm <- intersect(rf_genes, svm_genes)

cat("LASSO와 RF에서 겹치는 유전자:", common_lasso_rf, "\n")
cat("LASSO와 SVM-RFE에서 겹치는 유전자:", common_lasso_svm, "\n")
cat("RF와 SVM-RFE에서 겹치는 유전자:", common_rf_svm, "\n")

# 세 모델에서 겹치는 유전자
all_common_genes <- intersect(intersect(lasso_genes, rf_genes), svm_genes)
cat("세 모델에서 겹치는 유전자:", all_common_genes, "\n")














# ░▒▓ 패키지 로드 ▓▒░
library(pROC)
library(ggplot2)
library(RColorBrewer)

# ░▒▓ 1. ROC 객체 및 AUC 계산 ▓▒░
common_genes <- c("OASL", "CCDC142", "LTB", "TMEM229B", "BTN3A1", "RASAL3", "HP", "TMEM132B", "CD27", "AKAP12")

roc_list <- lapply(common_genes, function(g) {
  roc(test$group, test[[g]], quiet = TRUE, levels = rev(levels(factor(test$group))))
})
names(roc_list) <- common_genes
auc_vec <- sapply(roc_list, auc)

# ░▒▓ 2. AUC 기준 정렬 ▓▒░
ord <- order(auc_vec, decreasing = TRUE)
roc_list <- roc_list[ord]
auc_vec <- auc_vec[ord]
genes_sorted <- names(roc_list)

# ░▒▓ 3. 범례 라벨 및 색상 ▓▒░
legend_lab <- sprintf("%s (AUC = %.2f)", genes_sorted, auc_vec)
pal <- colorRampPalette(brewer.pal(9, "Set1"))(length(genes_sorted))

# ░▒▓ 4. ROC Curve 시각화 ▓▒░
ggroc(roc_list, aes = "color", legacy.axes = TRUE) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +  # 무작위 기준선
  scale_color_manual(values = pal,
                     breaks = genes_sorted,
                     labels = legend_lab,
                     name = "") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),               # 가운데 정렬 + 볼드
    legend.position = c(0.95, 0.05),                                     # 그래프 내부 오른쪽 아래
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.box.background = element_rect(color = "black", size = 0.3),
    legend.title = element_blank()
  ) +
  labs(title = "ROC Curve", x = "1 − Specificity", y = "Sensitivity")


# ░▒▓ AUC 바플롯: Top Genes by ROC AUC ▓▒░
library(ggplot2)

# AUC 정렬 정보로 데이터프레임 생성
auc_df <- data.frame(
  Gene = factor(names(auc_vec), levels = names(auc_vec)),  # 정렬된 순서 유지
  AUC = as.numeric(auc_vec)
)

ggplot(auc_df, aes(x = AUC, y = Gene)) +
  geom_bar(stat = "identity", fill = "lightsteelblue", color = "grey40", width = 0.7) +
  coord_cartesian(xlim = c(0.6, 1.00)) +  # 🔥 핵심 수정: 잘림 없이 시각적으로만 제한
  labs(title = "Top Genes by ROC AUC", x = "AUC", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )








# ---------------------
# 머신 러닝 + 공통 유전자 기반 ROC 분석
# --------------------
library(readr)
library(dplyr)
library(glmnet)
library(caret)
library(randomForest)
library(pROC)

# 1. 발현 데이터 불러오기 --------------------------------------------------
expr <- read_csv("GSE100927_gene_mapped_expression.csv", show_col_types = FALSE)
expr <- as.data.frame(expr)
rownames(expr) <- toupper(expr$symbol)
expr <- expr[, -1]  # symbol 컬럼 제거

# 2. 그룹 지정 -------------------------------------------------------------
group <- ifelse(grepl("Control", colnames(expr)), "NonDisease", "Disease")
group <- factor(group)
table(group)

# 3. CoDEG 유전자 리스트 불러오기 ------------------------------------------
genes3 <- toupper(scan("CoDEG3_UpUp_genes.txt", what = character()))
genes4 <- toupper(scan("CoDEG4_DownDown_genes.txt", what = character()))

# 4. 공통 분석 함수 정의 ---------------------------------------------------
run_common_roc_analysis <- function(gene_list, expr, group, label) {
  cat("\n🔍 [", label, "] 분석 시작...\n")
  
  # 발현 행렬 필터링 및 그룹 추가
  expr_sub <- expr[rownames(expr) %in% gene_list, , drop = FALSE]
  expr_sub <- t(expr_sub)
  expr_sub <- as.data.frame(expr_sub)
  expr_sub$group <- group
  
  # Train/Test 분할
  set.seed(123)
  train_idx <- createDataPartition(expr_sub$group, p = 0.7, list = FALSE)
  train <- expr_sub[train_idx, ]
  test  <- expr_sub[-train_idx, ]
  
  # ---------------- LASSO ------------------
  x_train <- train[, !colnames(train) %in% "group"] %>% mutate(across(everything(), as.numeric))
  x_train <- x_train[, apply(x_train, 2, var, na.rm = TRUE) != 0]
  x_matrix <- as.matrix(x_train)
  y_train <- ifelse(train$group == "Disease", 1, 0)
  x_df <- as.data.frame(x_matrix)
  x_df$group <- factor(ifelse(y_train == 1, "Disease", "NonDisease"), levels = c("NonDisease", "Disease"))
  
  ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
  lasso_model <- train(
    group ~ ., data = x_df,
    method = "glmnet",
    family = "binomial",
    metric = "ROC",
    trControl = ctrl,
    tuneGrid = expand.grid(alpha = 1, lambda = seq(0.0001, 1, length = 100))
  )
  coefs <- coef(lasso_model$finalModel, s = lasso_model$bestTune$lambda)
  lasso_genes <- rownames(coefs)[which(coefs != 0)]
  lasso_genes <- setdiff(lasso_genes, "(Intercept)")
  
  # ---------------- RF ------------------
  train_rf <- train
  train_rf$group <- factor(train_rf$group)
  x_rf <- train_rf[, !(colnames(train_rf) == "group")]
  y_rf <- train_rf$group
  rf_model <- randomForest(x = x_rf, y = y_rf, importance = TRUE)
  rf_importance <- importance(rf_model)
  rf_genes <- rownames(rf_importance)[order(rf_importance[, 1], decreasing = TRUE)][1:20]
  
  # ---------------- SVM-RFE ------------------
  ctrl_rfe <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
  svm_profile <- rfe(train[, -ncol(train)], train$group, sizes = c(5, 10, 20), rfeControl = ctrl_rfe)
  svm_genes <- predictors(svm_profile)
  
  # ---------------- 공통 유전자 추출 ------------------
  common_genes <- Reduce(intersect, list(lasso_genes, rf_genes, svm_genes))
  cat("✅ [", label, "] 세 모델 공통 유전자:", paste(common_genes, collapse=", "), "\n")
  
  # 두 모델에서 겹치는 유전자
  common_lasso_rf <- intersect(lasso_genes, rf_genes)
  common_lasso_svm <- intersect(lasso_genes, svm_genes)
  common_rf_svm <- intersect(rf_genes, svm_genes)
  
  cat("LASSO와 RF에서 겹치는 유전자:", common_lasso_rf, "\n")
  cat("LASSO와 SVM-RFE에서 겹치는 유전자:", common_lasso_svm, "\n")
  cat("RF와 SVM-RFE에서 겹치는 유전자:", common_rf_svm, "\n")
  
  # 세 모델에서 겹치는 유전자
  all_common_genes <- intersect(intersect(lasso_genes, rf_genes), svm_genes)
  cat("세 모델에서 겹치는 유전자:", all_common_genes, "\n")
  
  
  # ---------------- ROC 곡선 ------------------
  dir.create("ROC_Compare", showWarnings = FALSE)
  common_genes <- common_genes[common_genes %in% colnames(test)]
  if (length(common_genes) < 2) {
    message("❌ [", label, "] 공통 유전자가 2개 미만입니다. ROC 생략")
    return(NULL)
  }
  test_x <- test[, common_genes, drop = FALSE] %>% mutate(across(everything(), as.numeric))
  score <- rowMeans(test_x)
  group_bin <- ifelse(test$group == "Disease", 1, 0)
  
  roc_obj <- roc(group_bin, score, direction = ">")
  roc_path <- file.path("ROC_Compare", paste0(label, "_CommonML_ROC1.png"))
  cat("🖼️ ROC 이미지 저장 위치:", roc_path, "\n")
  png(roc_path, width = 600, height = 600)
  plot(roc_obj, col = "#e41a1c", lwd = 2, main = paste0(label, " Common ML ROC\nAUC = ", round(auc(roc_obj), 3)))
  abline(a = 0, b = 1, lty = 2, col = "gray")
  dev.off()
}

# ---------------- 실행 ------------------
run_common_roc_analysis(genes3, expr, group, label = "CoDEG3")
run_common_roc_analysis(genes4, expr, group, label = "CoDEG4")



# ░▒▓ 패키지 로드 ▓▒░
library(pROC)
library(ggplot2)
library(RColorBrewer)
library(patchwork)  # ROC + Barplot 같이 출력용

# ░▒▓ 1. 유전자 리스트 ▓▒░
common_genes <- c("BTN3A1", "CDKN2A", "CRHBP", "IFIT1", "LTB", "TMEM229B",
                  "IGFBP6", "MYBL1", "RGS5", "SCARA5", "SLPI", "UAP1")

# ░▒▓ 2. 그룹 지정 ▓▒░
group <- ifelse(grepl("Control", colnames(expr)), "NonDisease", "Disease")
group <- factor(group, levels = c("NonDisease", "Disease"))
group_bin <- ifelse(group == "Disease", 1, 0)

# ░▒▓ 3. ROC 객체 및 AUC 계산 ▓▒░
roc_list <- list()
auc_vec <- c()

for (gene in common_genes) {
  if (gene %in% rownames(expr)) {
    expr_vec <- as.numeric(expr[gene, ])
    roc_obj <- tryCatch(roc(group_bin, expr_vec, direction = ">"), error = function(e) NULL)
    if (!is.null(roc_obj)) {
      roc_list[[gene]] <- roc_obj
      auc_vec[gene] <- auc(roc_obj)
    }
  }
}

# ░▒▓ 4. AUC 기준 정렬 ▓▒░
ord <- order(auc_vec, decreasing = TRUE)
roc_list <- roc_list[ord]
auc_vec <- auc_vec[ord]
genes_sorted <- names(roc_list)

# ░▒▓ 5. 범례 라벨 및 색상 ▓▒░
legend_lab <- sprintf("%s (AUC = %.2f)", genes_sorted, auc_vec)
pal <- colorRampPalette(brewer.pal(9, "Set1"))(length(genes_sorted))

# ░▒▓ 6. ROC Curve 시각화 ▓▒░
p1 <- ggroc(roc_list, aes = "color", legacy.axes = TRUE) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
  scale_color_manual(values = pal,
                     breaks = genes_sorted,
                     labels = legend_lab,
                     name = "") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  ) +
  labs(title = "ROC Curve", x = "1 − Specificity", y = "Sensitivity")

# ░▒▓ 7. AUC 바플롯 ▓▒░
auc_df <- data.frame(
  Gene = factor(names(auc_vec), levels = names(auc_vec)),
  AUC = as.numeric(auc_vec)
)

p2 <- ggplot(auc_df, aes(x = AUC, y = Gene)) +
  geom_bar(stat = "identity", fill = "lightsteelblue", color = "grey40", width = 0.7) +
  coord_cartesian(xlim = c(0.6, 1.00)) +
  labs(title = "Top Genes by ROC AUC", x = "AUC", y = "") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# ░▒▓ 8. 한 화면에 ROC + AUC plot ▓▒░
p1
p2







# 0. 필수 패키지 로드
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 1. CoDEG3 데이터 읽기
deg3 <- read_csv("CoDEG3_Annotated.csv", show_col_types = FALSE)

# 2. UpUp/DownDown 유전자만 추출
genes3 <- deg3 %>%
  filter(direction %in% c("UpUp","DownDown")) %>%
  pull(gene) %>%
  toupper() %>%
  unique()

# 3. SYMBOL → ENTREZID 변환
ent3 <- bitr(genes3,
             fromType = "SYMBOL",
             toType   = "ENTREZID",
             OrgDb    = org.Hs.eg.db) %>%
  distinct(ENTREZID) %>%
  pull(ENTREZID)

# 4. GO BP enrichment
ego3 <- enrichGO(gene         = ent3,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 pAdjustMethod= "BH",
                 pvalueCutoff = 0.05,
                 readable     = TRUE)

# 5. 플롯에서 사용된 주요 BP 용어들
key_terms3 <- c(
  "positive regulation of cell adhesion",
  "leukocyte cell-cell adhesion",
  "regulation of T cell activation",
  "regulation of leukocyte cell-cell adhesion",
  "leukocyte proliferation",
  "positive regulation of leukocyte activation",
  "mononuclear cell proliferation",
  "lymphocyte proliferation",
  "immune response-regulating cell surface receptor signaling pathway",
  "T cell proliferation"
)

# 6. 해당 용어에 매핑된 GO 결과 필터
sel3 <- ego3@result %>%
  filter(Description %in% key_terms3)

# 7. 각 term에 속한 유전자 추출
genes_for_ppi3 <- sel3$geneID %>%     # "GENE1/GENE2/..." 형식
  strsplit("/") %>%
  unlist() %>%
  unique()

# 8. 결과 출력
cat("PPI 네트워크용 CoDEG3 유전자 (GO BP 기반):\n")
print(genes_for_ppi3)

# -------------------------------
# CoDEG3 genes_for_ppi3 기반 PPI 유사 네트워크 분석
# -------------------------------

# 0. 필수 패키지 로드
library(dplyr)
library(purrr)
library(igraph)

# 1. genes_for_ppi3: 이전 단계에서 얻은 유전자 벡터
#    ex) genes_for_ppi3 <- c("CD74","CTSS","HLA-DRA",...)

# 2. 발현행렬 로드 및 전처리
expr_raw <- read.csv(
  "GSE100927_gene_mapped_expression.csv",
  row.names   = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expr_mat <- t(expr_raw)                    # 행=샘플, 열=유전자
colnames(expr_mat) <- toupper(colnames(expr_mat))

# 3. 사용할 유전자 매칭
genes3 <- toupper(genes_for_ppi3)
matched <- intersect(genes3, colnames(expr_mat))
if (length(matched) < 3) stop("매칭 유전자 수 부족: 확인 필요")

expr_sub <- expr_mat[, matched]

# 4. 유전자 간 상관행렬 계산
cor_mtx <- cor(expr_sub, use = "pairwise.complete.obs")

# 5. long format 변환 및 자기자신 제거
cor_df <- as.data.frame(as.table(cor_mtx))
colnames(cor_df) <- c("from","to","cor")
cor_df <- cor_df %>% filter(from != to)

# 6. 각 유전자별 상위 3개 엣지 추출
top_edges <- map_dfr(matched, function(g) {
  cor_df %>%
    filter(from == g | to == g) %>%
    arrange(desc(abs(cor))) %>%
    slice_head(n = 3)
})

# 7. igraph 객체 생성 및 허브 유전자 추출
g3 <- graph_from_data_frame(top_edges, directed = FALSE)

hub_score3   <- degree(g3, mode = "all")
top_hubs3    <- names(sort(hub_score3, decreasing = TRUE))[1:10]

# 8. 시각화 속성 설정
V(g3)$color      <- ifelse(V(g3)$name %in% top_hubs3, "tomato", "skyblue")
V(g3)$size       <- ifelse(V(g3)$name %in% top_hubs3, 25, 12)
V(g3)$label.cex  <- ifelse(V(g3)$name %in% top_hubs3, 1.6, 1.0)
V(g3)$label.font <- ifelse(V(g3)$name %in% top_hubs3, 2, 1)
V(g3)$label.color<- "black"

E(g3)$width <- abs(E(g3)$cor) * 4
E(g3)$color <- "gray40"

# 9. 레이아웃 생성 및 플롯
set.seed(42)
lay3 <- layout_with_fr(g3)

plot(g3,
     layout              = lay3,
     vertex.label.family = "sans",
     edge.curved         = 0.1,
     main                = "CoDEG3 GO-BP / PPI Network")

# 1. PNG로 저장
png("CoDEG3_PPI_Network.png",
    width  = 10, height = 10,
    units  = "in", res    = 300)
plot(g3,
     layout              = lay3,
     vertex.label.family = "sans",
     edge.curved         = 0.1,
     main                = "CoDEG3 GO-BP / PPI Network")
dev.off()







# 0. 필수 패키지 로드
library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# 1. CoDEG4 데이터 읽기
deg4 <- read_csv("CoDEG4_Annotated.csv", show_col_types = FALSE)

# 2. UpUp/DownDown 유전자만 추출
genes4 <- deg4 %>%
  filter(direction %in% c("UpUp","DownDown")) %>%
  pull(gene) %>%
  toupper() %>%
  unique()

# 3. SYMBOL → ENTREZID 변환
ent4 <- bitr(genes4,
             fromType = "SYMBOL",
             toType   = "ENTREZID",
             OrgDb    = org.Hs.eg.db) %>%
  distinct(ENTREZID) %>%
  pull(ENTREZID)

# 4. GO BP enrichment
ego4 <- enrichGO(gene         = ent4,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 pAdjustMethod= "BH",
                 pvalueCutoff = 0.05,
                 readable     = TRUE)

# 5. 주요 BP term 정의 (플롯에서 확인된 것들)
key_terms <- c(
  "positive regulation of cell activation",
  "positive regulation of leukocyte activation",
  "antigen processing and presentation of peptide antigen via MHC class II",
  "collagen metabolic process"
)

# 6. 위 term에 해당하는 GO 결과만 필터
sel_bp <- ego4@result %>%
  filter(Description %in% key_terms)

# 7. 각 term에 속한 유전자 추출
genes_for_ppi <- sel_bp$geneID %>%     # "CD74/CTSS/HLA-DRA/..." 형식
  strsplit("/") %>%
  unlist() %>%
  unique()

# 8. 결과 확인
cat("PPI 네트워크용 CoDEG4 유전자 (GO BP 기반):\n")
print(genes_for_ppi)

# → 이 `genes_for_ppi` 벡터를 igraph/plot 함수의 node 리스트로 사용하세요.

# -------------------------------
# CoDEG4 genes_for_ppi 기반 PPI 유사 네트워크 분석
# -------------------------------

# 0. 필수 패키지 로드
library(dplyr)
library(purrr)
library(igraph)

# 1. genes_for_ppi: 이전 단계에서 얻은 CoDEG4용 유전자 벡터
#    ex) genes_for_ppi <- c("COL1A1","COL3A1","CD74",...)

# 2. 발현행렬 로드 및 전처리
expr_raw <- read.csv(
  "GSE100927_gene_mapped_expression.csv",
  row.names   = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expr_mat <- t(expr_raw)                     # 행=샘플, 열=유전자
colnames(expr_mat) <- toupper(colnames(expr_mat))

# 3. 사용할 유전자 매칭
genes4 <- toupper(genes_for_ppi)
matched4 <- intersect(genes4, colnames(expr_mat))
if (length(matched4) < 3) stop("매칭 유전자 수가 너무 적습니다. 확인 필요.")

expr_sub4 <- expr_mat[, matched4]

# 4. 유전자 간 상관행렬 계산
cor_mtx4 <- cor(expr_sub4, use = "pairwise.complete.obs")

# 5. long format 변환 및 자기자신 제거
cor_df4 <- as.data.frame(as.table(cor_mtx4))
colnames(cor_df4) <- c("from", "to", "cor")
cor_df4 <- cor_df4 %>% filter(from != to)

# 6. 각 유전자별 상위 3개 엣지 추출
top_edges4 <- map_dfr(matched4, function(g) {
  cor_df4 %>%
    filter(from == g | to == g) %>%
    arrange(desc(abs(cor))) %>%
    slice_head(n = 3)
})

# 7. igraph 객체 생성 및 허브 유전자 추출
g4          <- graph_from_data_frame(top_edges4, directed = FALSE)
hub_scores4 <- degree(g4, mode = "all")
top_hubs4   <- names(sort(hub_scores4, decreasing = TRUE))[1:10]

# 8. 시각화 속성 설정
V(g4)$color      <- ifelse(V(g4)$name %in% top_hubs4, "tomato", "skyblue")
V(g4)$size       <- ifelse(V(g4)$name %in% top_hubs4, 25, 12)
V(g4)$label.cex  <- ifelse(V(g4)$name %in% top_hubs4, 1.6, 1.0)
V(g4)$label.font <- ifelse(V(g4)$name %in% top_hubs4, 2, 1)
V(g4)$label.color<- "black"
E(g4)$width      <- abs(E(g4)$cor) * 4
E(g4)$color      <- "gray40"

# 9. 레이아웃 생성 및 네트워크 플롯
set.seed(99)
lay4 <- layout_with_fr(g4)

plot(g4,
     layout              = lay4,
     vertex.label.family = "sans",
     edge.curved         = 0.1,
     main                = "CoDEG4 GO-BP / PPI Network")



# 9. 파일로 저장 (PNG)
png("CoDEG4_PPI_Network.png",
    width  = 10, height = 10,
    units  = "in", res    = 300)
plot(g4,
     layout              = lay4,
     vertex.label.family = "sans",
     edge.curved         = 0.1,
     main                = "CoDEG4 GO-BP / PPI Network")
dev.off()


# CoDEG3 네트워크(g3)에서 허브 유전자 도출
# 1. degree 계산
hub_score3 <- degree(g3, mode = "all")

# 2. 내림차순 정렬
hub_score3_sorted <- sort(hub_score3, decreasing = TRUE)

# 3. 상위 10개 허브 유전자 추출
top_hubs3 <- head(hub_score3_sorted, 11)

# 4. 결과 출력
print(top_hubs3)

# CoDEG4 네트워크(g4)에서 허브 유전자 도출
# 1. degree 계산
hub_score4 <- degree(g4, mode = "all")

# 2. 내림차순 정렬
hub_score4_sorted <- sort(hub_score4, decreasing = TRUE)

# 3. 상위 10개 허브 유전자 추출
top_hubs4 <- head(hub_score4_sorted, 11)

# 4. 결과 출력
print(top_hubs4)



# 패키지 로드
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# 1. Seurat 객체는 'seurat'로 가정
Idents(seurat) <- "cell_type_singler"

# 2. 허브 유전자 리스트
hub_list <- list(
  CoDEG3 = c("SKAP1","TYROBP","CD14","CCL5","LIPA","FBLN1",
             "CD3G","HLA-DMA","CD5","RASAL3","THEMIS"),
  CoDEG4 = c("HLA-DMA","HLA-DRB1","CTSB","TYROBP","COL5A1",
             "PCOLCE","HLA-DMB","CCL3","CD74","IL7R","CD83")
)

# 3. OASL 그룹 변수 확인
if (!"OASL_group" %in% colnames(seurat@meta.data)) {
  stop("메타데이터에 'OASL_group' 컬럼이 없습니다.")
}

# 4. 면역세포 타입 정의
immune_types <- c("CD4+ T-cells","CD8+ T-cells","NK cells",
                  "B-cells","DC","Monocytes","Macrophages")

for (cond in names(hub_list)) {
  hub_genes <- hub_list[[cond]]
  
  # A) UMAP
  p_umap <- FeaturePlot(
    seurat,
    features = hub_genes,
    reduction = "umap",
    cols      = c("gray90","red"),
    pt.size   = 0.3
  ) + ggtitle(paste(cond, "UMAP")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(cond, "_UMAP.png"), p_umap, width=10, height=6, dpi=300)
  
  # B) DotPlot
  p_dot <- DotPlot(
    seurat,
    features = hub_genes,
    group.by = "cell_type_singler"
  ) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste(cond, "DotPlot"))
  ggsave(paste0(cond, "_DotPlot.png"), p_dot, width=8, height=6, dpi=300)
  
  # C) RidgePlot
  Idents(seurat) <- "OASL_group"
  p_ridge <- RidgePlot(
    seurat,
    features = hub_genes,
    group.by = "OASL_group",
    ncol     = 2
  ) + ggtitle(paste(cond, "RidgePlot (OASL High vs Low)")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(cond, "_RidgePlot.png"), p_ridge, width=10, height=8, dpi=300)
  
  # D) 면역세포 ViolinPlot
  seu_immune <- subset(seurat, subset = cell_type_singler %in% immune_types)
  Idents(seu_immune) <- "OASL_group"
  p_violin_list <- VlnPlot(
    seu_immune,
    features = hub_genes,
    pt.size  = 0,
    group.by = "OASL_group",
    split.by = "cell_type_singler",
    combine  = FALSE
  )
  p_violin <- wrap_plots(p_violin_list, ncol = 2) +
    plot_annotation(title = paste(cond, "Immune Cells Violin (OASL High vs Low)")) &
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0(cond, "_Immune_Violin.png"), p_violin, width=12, height=10, dpi=300)
  
  # 다음 조건을 위해 identity 초기화
  Idents(seurat) <- "cell_type_singler"
}


##### 로지스틱 회귀 모델 #####
# 1. Signature score 계산
# CoDEG3, CoDEG4 유전자 벡터 (대문자)
genes3 <- toupper(c("BTN3A1","CDKN2A","CRHBP","IFIT1","IGFBP6","LTB"))  # 앞서 뽑은 CoDEG3 허브
genes4 <- toupper(c("MYBL1","RGS5","SCARA5","SLPI","TMEM229B","UAP1"))  # CoDEG4 허브

# 매트릭스에서 매칭
mat3 <- expr_mat[, intersect(genes3, colnames(expr_mat))]
mat4 <- expr_mat[, intersect(genes4, colnames(expr_mat))]

# 평균 발현으로 signature score 생성
score3 <- rowMeans(mat3)
score4 <- rowMeans(mat4)

# 데이터프레임 결합
df <- data.frame(
  group  = factor(group, levels=c("NonDisease","Disease")),
  score3 = score3,
  score4 = score4
)

# 2. 단일 서명 logistic 회귀
# CoDEG3 score
fit3 <- glm(group ~ score3, data=df, family="binomial")
summary(fit3)            # 회귀계수, p-value 확인
library(pROC)
roc3 <- roc(df$group, predict(fit3, type="response"))
auc3 <- auc(roc3); print(auc3)

# CoDEG4 score
fit4 <- glm(group ~ score4, data=df, family="binomial")
summary(fit4)
roc4 <- roc(df$group, predict(fit4, type="response"))
auc4 <- auc(roc4); print(auc4)


# 3. 서명 결합 모델
fit34 <- glm(group ~ score3 + score4, data=df, family="binomial")
summary(fit34)
roc34 <- roc(df$group, predict(fit34, type="response"))
auc34 <- auc(roc34); print(auc34)


# 4. 다중 유전자 glm (페널티 없이)
# CoDEG3 유전자 개별 회귀
df3g <- data.frame(group=df$group, mat3)
fit3g <- glm(group ~ ., data=df3g, family="binomial")
summary(fit3g)  # 각 유전자 β, 유의성 확인

# CoDEG4 동일
df4g <- data.frame(group=df$group, mat4)
fit4g <- glm(group ~ ., data=df4g, family="binomial")
summary(fit4g)


######## 시각화 #########
######## 1. ROC 곡선 비교 플롯 및 저장 ########
library(pROC); library(ggplot2)

roc3  <- roc(df$group, predict(fit3,  type="response"))
roc4  <- roc(df$group, predict(fit4,  type="response"))
roc34 <- roc(df$group, predict(fit34, type="response"))

roc_df <- rbind(
  data.frame( sens=rev(roc3$sensitivities),  spec=rev(roc3$specificities),  model="CoDEG3"),
  data.frame( sens=rev(roc4$sensitivities),  spec=rev(roc4$specificities),  model="CoDEG4"),
  data.frame( sens=rev(roc34$sensitivities), spec=rev(roc34$specificities), model="Combined")
)

p_roc <- ggplot(roc_df, aes(x=1-spec, y=sens, color=model)) +
  geom_line(size=1.2) +
  geom_abline(linetype="dashed") +
  labs(x="False Positive Rate", y="True Positive Rate",
       title="ROC Curve Comparison",
       subtitle=paste0("AUC: CoDEG3=", round(auc(roc3),3),
                       ", CoDEG4=", round(auc(roc4),3),
                       ", Combined=", round(auc(roc34),3))) +
  scale_color_manual(values=c("CoDEG3"="firebrick","CoDEG4"="darkgreen","Combined"="steelblue")) +
  theme_classic(base_size=14) +
  theme(
    plot.title     = element_text(face="bold"),
    plot.subtitle  = element_text(size=12),
    legend.position="bottom",
    panel.border   = element_rect(color="black", fill=NA),
    axis.ticks     = element_line(color="black")
  )

# 저장
ggsave("Figure_ROC_Comparison.png", p_roc, width=8, height=6, dpi=300)


######## 2. 회귀계수 포레스트 플롯 및 저장 ########
library(broom); library(ggplot2)

coef_df <- tidy(fit34) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    odds_ratio = exp(estimate),
    lower      = exp(estimate - 1.96*std.error),
    upper      = exp(estimate + 1.96*std.error)
  )

p_forest <- ggplot(coef_df, aes(x=term, y=odds_ratio)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2) +
  geom_hline(yintercept=1, linetype="dashed", color="gray40") +
  coord_flip() +
  labs(x=NULL, y="Odds Ratio (95% CI)",
       title="Coefficients of Combined Logistic Model") +
  theme_classic(base_size=14) +
  theme(
    plot.title   = element_text(face="bold"),
    panel.border = element_rect(color="black", fill=NA),
    axis.ticks   = element_line(color="black")
  )

# 저장
ggsave("Figure_Forest_Coefficients.png", p_forest, width=6, height=4, dpi=300)


######## 3. Signature Score 분포 박스/바이올린 플롯 및 저장 ########
library(reshape2); library(ggplot2)

df_m <- melt(df, id.vars="group", measure.vars=c("score3","score4"))

p_score <- ggplot(df_m, aes(x=group, y=value, fill=group)) +
  geom_violin(alpha=0.6) +
  geom_boxplot(width=0.2, outlier.shape=NA) +
  facet_wrap(~variable, scales="free_y") +
  labs(x="", y="Signature Score",
       title="Signature Score Distribution by Group") +
  scale_fill_manual(values=c("NonDisease"="#E69F00","Disease"="#56B4E9")) +
  theme_classic(base_size=14) +
  theme(
    plot.title    = element_text(face="bold"),
    legend.position="none",
    panel.border  = element_rect(color="black", fill=NA),
    axis.ticks    = element_line(color="black")
  )

# 저장
ggsave("Figure_Score_Distribution.png", p_score, width=8, height=4, dpi=300)


######## 4. Calibration Curve 및 저장 ########
library(dplyr); library(ggplot2)

df <- df %>% mutate(pred34 = predict(fit34, type = "response"))
cal_df <- df %>%
  mutate(bin = ntile(pred34, 10)) %>%
  group_by(bin) %>%
  summarize(pred = mean(pred34), obs = mean(group=="Disease"))

p_cal <- ggplot(cal_df, aes(x=pred, y=obs)) +
  geom_line(size=1.2, color="steelblue") +
  geom_point(size=3, color="steelblue") +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray40") +
  labs(x="Predicted Probability", y="Observed Proportion",
       title="Calibration Curve (Combined Model)") +
  theme_classic(base_size=14) +
  theme(
    plot.title   = element_text(face="bold", hjust=0.5),
    panel.border = element_rect(color="black", fill=NA),
    axis.ticks   = element_line(color="black")
  )

# 저장
ggsave("Figure_Calibration_Curve.png", p_cal, width=6, height=6, dpi=300)






######## 1. 패키지 로딩 ########
library(pROC)
library(ggplot2)
library(dplyr)

######## 2. OASL 발현값 추가 (expr_mat: 샘플 x 유전자 구조) ########
df$OASL <- as.numeric(expr_mat[, "OASL"])

######## 3. 로지스틱 회귀 모델 생성 ########
fit_oasl <- glm(group ~ OASL, data = df, family = "binomial")             # OASL 단독
fit3     <- glm(group ~ score3, data = df, family = "binomial")           # CoDEG3 단독
fit4     <- glm(group ~ score4, data = df, family = "binomial")           # CoDEG4 단독
fit34    <- glm(group ~ score3 + score4, data = df, family = "binomial")  # 3+4 결합
fit_all  <- glm(group ~ score3 + score4 + OASL, data = df, family = "binomial") # 최종 결합

######## 4. ROC 객체 생성 ########
roc_oasl <- roc(df$group, predict(fit_oasl, type = "response"))
roc3     <- roc(df$group, predict(fit3, type = "response"))
roc4     <- roc(df$group, predict(fit4, type = "response"))
roc34    <- roc(df$group, predict(fit34, type = "response"))
roc_all  <- roc(df$group, predict(fit_all, type = "response"))

######## 5. AUC 포함 모델 이름 생성 ########
label_oasl <- paste0("OASL (AUC=", round(auc(roc_oasl), 3), ")")
label_3    <- paste0("CoDEG3 (AUC=", round(auc(roc3), 3), ")")
label_4    <- paste0("CoDEG4 (AUC=", round(auc(roc4), 3), ")")
label_34   <- paste0("Combined_34 (AUC=", round(auc(roc34), 3), ")")
label_all  <- paste0("Combined_All (AUC=", round(auc(roc_all), 3), ")")

######## 6. ROC 데이터프레임 생성 ########
roc_df <- rbind(
  data.frame(sens = rev(roc_oasl$sensitivities), spec = rev(roc_oasl$specificities), model = label_oasl),
  data.frame(sens = rev(roc3$sensitivities),     spec = rev(roc3$specificities),     model = label_3),
  data.frame(sens = rev(roc4$sensitivities),     spec = rev(roc4$specificities),     model = label_4),
  data.frame(sens = rev(roc34$sensitivities),    spec = rev(roc34$specificities),    model = label_34),
  data.frame(sens = rev(roc_all$sensitivities),  spec = rev(roc_all$specificities),  model = label_all)
)

######## 7. 색상 매핑 정의 ########
color_mapping <- setNames(
  c("orange", "firebrick", "darkgreen", "steelblue", "purple"),
  c(label_oasl, label_3, label_4, label_34, label_all)
)

######## 8. ROC 곡선 시각화 ########
p_roc <- ggplot(roc_df, aes(x = 1 - spec, y = sens, color = model)) +
  geom_line(size = 1.2) +
  geom_abline(linetype = "dashed") +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = "ROC Curve Comparison") +
  scale_color_manual(
    name = "Model (AUC)",
    values = color_mapping
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title     = element_text(face = "bold", size = 18),   # 제목 크게
    legend.position = "bottom",
    legend.title    = element_text(face = "bold", size = 14),  # 범례 제목 크게
    legend.text     = element_text(size = 10),                 # 범례 항목 작게
    panel.border    = element_rect(color = "black", fill = NA),
    axis.ticks      = element_line(color = "black")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))  # 범례 두 줄 정렬

######## 9. 저장 ########
ggsave("Figure_ROC_Comparison_All_with_Legend.png", p_roc, width = 8, height = 6, dpi = 300)

#  3. Signature Score + OASL 분포 (박스+바이올린)
######## 3. Signature Score + OASL 분포 플롯 ########
library(reshape2)
library(ggplot2)

# score3, score4, OASL 포함 melt
df_m <- melt(df, id.vars = "group", measure.vars = c("score3", "score4", "OASL"))

# 플롯
p_score <- ggplot(df_m, aes(x = group, y = value, fill = group)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", size = 0.4) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = NULL, y = "Signature Value",
       title = "Distribution of Signature Scores by Group") +
  scale_fill_manual(values = c("NonDisease" = "#E69F00", "Disease" = "#56B4E9")) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    strip.text    = element_text(face = "bold"),
    legend.position = "none",
    panel.border  = element_rect(color = "black", fill = NA),
    axis.ticks    = element_line(color = "black")
  )

# 저장
ggsave("Figure_Score_Distribution1.png", p_score, width = 9, height = 4.5, dpi = 300)


# 4. Calibration Curve (최종 결합 모델)
######## 4. Calibration Curve (OASL + score3 + score4) ########
library(dplyr)
library(ggplot2)

# 예측값 추가
df <- df %>% mutate(pred_all = predict(fit_all, type = "response"))

# 10등분 binning
cal_df <- df %>%
  mutate(bin = ntile(pred_all, 10)) %>%
  group_by(bin) %>%
  summarize(pred = mean(pred_all), obs = mean(group == "Disease"))

# 플롯
p_cal <- ggplot(cal_df, aes(x = pred, y = obs)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(x = "Predicted Probability", y = "Observed Proportion",
       title = "Calibration Curve (OASL + CoDEG3 + CoDEG4)") +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks   = element_line(color = "black")
  )

# 저장
ggsave("Figure_Calibration_Curve1.png", p_cal, width = 6, height = 6, dpi = 300)


