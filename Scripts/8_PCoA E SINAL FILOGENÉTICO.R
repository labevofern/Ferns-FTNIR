# =============================================================================
# PCoA E SINAL FILOGENÉTICO
# =============================================================================
# Article: ........................... 
# Authors: Marise H.V. de Oliveira, Niksoney A. Mendonça, D. Victor S. Silva, 
#          Julia Z. Alves, Keven dos S. Lima, Carolina do Valle, 
#          Juliana Aljahara and Thaís E. Almeida
# Script author: Niksoney Azevedo Mendonça
# E-mail: niksoneyazevedo2017@gmail.com
# =============================================================================
# PACOTES PARA NECESSÁRIOS
# =============================================================================
install.packages(c("ape","phytools","dplyr","ggplot2","ggrepel","phylosignal","phylobase", "patchword"))

library(ape)           # análises e manipulação de árvores filogenéticas
library(phytools)      # métodos e visualizações filogenéticas
library(dplyr)         # manipulação de dados
library(ggplot2)       # visualização de dados
library(ggrepel)       # labels sem sobreposição
library(phylosignal)   # cálculo de sinal filogenético
library(phylobase)     # estrutura de dados filogenéticos
library(patchwork)     # organização de múltiplos gráficos
library(viridis)       # escalas de cores amigáveis
# =============================================================================
# Definir diretório
# =============================================================================
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/Sinal filogenético")

# =============================================================================
# Ler arvore e renomear
# =============================================================================
tree_file <- "ftol_subtree_especies_selecionadas.tree"
tree <- read.tree(tree_file)

renomear_Species <- data.frame(
  old = c("DICFLE","LYGVEN","ANEHIR","CYAMIC","DIPCRI",
          "BLEOCC","MENSER","CYCMEN","PLEAST","MICGEM",
          "NIPCRA","SERTRI","ADIPUL"),
  new = c("Dicranopteris flexuosa","Lygodium volubile","Anemia hirta","Cyathea microdonta","Diplazium cristatum",
          "Blechnum occidentale","Meniscium serratum","Cyclodium meniscioides","Pleoeltis astrolepis","Microgramma geminata",
          "Niphidium crassifolium","Serpocaulon triseriale","Adiantum pulverulentum")
)

for(i in 1:nrow(renomear_Species)){
  tree$tip.label[tree$tip.label == renomear_Species$old[i]] <- renomear_Species$new[i]
}

# =============================================================================
# Ler dados dos espectros médios
# =============================================================================
abaxial_spectrum <- read.csv("Matriz_abaxial_media_todas os espectros.csv", stringsAsFactors = FALSE)
adaxial_spectrum <- read.csv("Matriz_adaxial_media_todas os espectros.csv", stringsAsFactors = FALSE)
combined <- read.csv("Matriz_Combined_media_todas os espectros.csv", stringsAsFactors = FALSE)

# =============================================================================
# Calcular matriz de distância euclidiana
# =============================================================================
calcular_distancia <- function(spectrum_data){
  rownames(spectrum_data) <- spectrum_data[,1]
  spectrum_data <- spectrum_data[,-1]
  dist(spectrum_data, method = "euclidean")
}

dist_abaxial <- calcular_distancia(abaxial_spectrum)
dist_adaxial <- calcular_distancia(adaxial_spectrum)
dist_combined <- calcular_distancia(combined)

# =============================================================================
# Fazer PcoA
# =============================================================================
fazer_pcoa <- function(dist_matrix, k = 2, tree = NULL){
  pcoa_result <- cmdscale(dist_matrix, k = k, eig = TRUE)
  coords <- as.data.frame(pcoa_result$points)
  colnames(coords) <- paste0("PCoA", 1:k)
  if(!is.null(tree)) rownames(coords) <- tree$tip.label
  coords
}

pcoa_abaxial <- fazer_pcoa(dist_abaxial, tree = tree)
pcoa_adaxial <- fazer_pcoa(dist_adaxial, tree = tree)
pcoa_combined <- fazer_pcoa(dist_combined, tree = tree)

# =============================================================================
# Calcular sinal filogenético K
# =============================================================================
calcular_K <- function(coords, tree){
  K_list <- list()
  for(i in 1:ncol(coords)){
    trait <- coords[,i]
    names(trait) <- rownames(coords)
    K_list[[colnames(coords)[i]]] <- phylosig(tree, trait, method = "K", test = TRUE)
  }
  K_list
}

K_abaxial <- calcular_K(pcoa_abaxial, tree)
K_adaxial <- calcular_K(pcoa_adaxial, tree)
K_combined <- calcular_K(pcoa_combined, tree)

# =============================================================================
#  Tabela com valores de K
# =============================================================================
extrair_K_info <- function(K_list, spectrum_name){
  df <- data.frame(
    PCoA = names(K_list),
    K_value = sapply(K_list, function(x) x$K),
    P_value = sapply(K_list, function(x) x$P),
    Spectrum = spectrum_name
  )
  df
}

K_table_abaxial <- extrair_K_info(K_abaxial, "Abaxial")
K_table_adaxial <- extrair_K_info(K_adaxial, "Adaxial")
K_table_combined <- extrair_K_info(K_combined, "Combined")

K_table_all <- rbind(K_table_abaxial, K_table_adaxial, K_table_combined)
write.csv(K_table_all, "Sinal_filogenetico_K_tabela.csv", row.names = FALSE)

# =============================================================================
# Tabela por espécie com valores de PCoA
# =============================================================================
criar_tabela_por_especie <- function(pcoa_coords, spectrum_name){
  df <- pcoa_coords
  df$Species <- rownames(df)
  df$Spectrum <- spectrum_name
  df
}

pcoa_table_abaxial <- criar_tabela_por_especie(pcoa_abaxial, "Abaxial")
pcoa_table_adaxial <- criar_tabela_por_especie(pcoa_adaxial, "Adaxial")
pcoa_table_combined <- criar_tabela_por_especie(pcoa_combined, "Combined")

pcoa_table_all <- rbind(pcoa_table_abaxial, pcoa_table_adaxial, pcoa_table_combined)
write.csv(pcoa_table_all, "PCoA_por_especie.csv", row.names = FALSE)

# =============================================================================
# Visualização Filogenética - heatmap PCoA1
# =============================================================================

heatmap_mat <- data.frame(
  Abaxial = pcoa_abaxial$PCoA1,
  Adaxial = pcoa_adaxial$PCoA1,
  Combined = pcoa_combined$PCoA1
)
rownames(heatmap_mat) <- rownames(pcoa_abaxial)  # garantir que bate com tree$tip.label

png("heatmap_filogenetico_PCoA1.png",
    width = 14, height = 8, units = "in", res = 300)

par(mar = c(8, 6, 8, 4) + 0.1)  # Aumentei as margens para caber textos maiores

# Plot principal com textos aumentados
phylo.heatmap(tree, heatmap_mat,
              standardize = TRUE,
              fsize = 1.4,        # Aumentei o tamanho da fonte da árvore
              length = 0.8)

# Calcular posições corretas para os nomes
n_cols <- ncol(heatmap_mat)
posicoes <- seq(0.5, n_cols - 0.5, by = 1)

# Adicionar nomes retos acima do heatmap com fonte maior
mtext(text = colnames(heatmap_mat), 
      side = 3, 
      at = posicoes,
      line = 2,                   # Aumentei a distância
      cex = 1.3,                  # Aumentei o tamanho da fonte
      las = 1)

# Se quiser também aumentar os nomes das linhas (espécies)
# Adicionar esta linha se a função phylo.heatmap permitir controlar labels laterais
mtext(text = rownames(heatmap_mat), 
      side = 4, 
      at = 1:nrow(heatmap_mat),
      line = 1, 
      cex = 1.2,                  # Tamanho para nomes laterais
      las = 1)

dev.off()

# =============================================================================
# Visualização Filogenética - heatmap PCoA2
# =============================================================================

heatmap_mat <- data.frame(
  Abaxial = pcoa_abaxial$PCoA2,
  Adaxial = pcoa_adaxial$PCoA2,
  Combined = pcoa_combined$PCoA2
)
rownames(heatmap_mat) <- rownames(pcoa_abaxial)  # garantir que bate com tree$tip.label

png("heatmap_filogenetico_PCoA2.png",
    width = 14, height = 8, units = "in", res = 300)

par(mar = c(8, 6, 8, 4) + 0.1)  # Aumentei as margens para caber textos maiores

# Plot principal com textos aumentados
phylo.heatmap(tree, heatmap_mat,
              standardize = TRUE,
              fsize = 1.4,        # Aumentei o tamanho da fonte da árvore
              length = 0.8)

# Calcular posições corretas para os nomes
n_cols <- ncol(heatmap_mat)
posicoes <- seq(0.5, n_cols - 0.5, by = 1)

# Adicionar nomes retos acima do heatmap com fonte maior
mtext(text = colnames(heatmap_mat), 
      side = 3, 
      at = posicoes,
      line = 2,                   # Aumentei a distância
      cex = 1.3,                  # Aumentei o tamanho da fonte
      las = 1)

# Se quiser também aumentar os nomes das linhas (espécies)
# Adicionar esta linha se a função phylo.heatmap permitir controlar labels laterais
mtext(text = rownames(heatmap_mat), 
      side = 4, 
      at = 1:nrow(heatmap_mat),
      line = 1, 
      cex = 1.2,                  # Tamanho para nomes laterais
      las = 1)

dev.off()

# =============================================================================
# Gráficos da PCoA e Scree
# =============================================================================
renomear_Species <- data.frame(
  new = c("DICFLE","LYGVEN","ANEHIR","CYAMIC","DIPCRI",
          "BLEOCC","MENSER","CYCMEN","PLEAST","MICGEM",
          "NIPCRA","SERTRI","ADIPUL"),
  old = c("Dicranopteris flexuosa","Lygodium volubile","Anemia hirta","Cyathea microdonta","Diplazium cristatum",
          "Blechnum occidentale","Meniscium serratum","Cyclodium meniscioides","Pleopeltis astrolepis","Microgramma geminata",
          "Niphidium crassifolium","Serpocaulon triseriale","Adiantum pulverulentum"),
  stringsAsFactors = FALSE
)

# utilidade: formatar e mapear labels
format_species <- function(x){
  x2 <- gsub("_", " ", x)
  trimws(x2)
}
map_to_new_label <- function(label, ren_df){
  idx <- match(label, ren_df$old)
  if(!is.na(idx)) return(ren_df$new[idx])
  idx2 <- match(tolower(label), tolower(ren_df$old))
  if(!is.na(idx2)) return(ren_df$new[idx2])
  label
}

# ---------------------------
# PCoA + utilidades
# ---------------------------
fazer_pcoa <- function(dist_matrix, k = 10, tree = NULL, ren_df = NULL){
  pcoa_result <- cmdscale(dist_matrix, k = k, eig = TRUE)
  coords <- as.data.frame(pcoa_result$points)
  colnames(coords) <- paste0("PCoA", 1:k)
  if(!is.null(tree)) rownames(coords) <- tree$tip.label
  coords$species <- rownames(coords)
  coords$label_orig <- format_species(coords$species)
  if(!is.null(ren_df)){
    coords$label <- vapply(coords$label_orig, map_to_new_label, FUN.VALUE = "", ren_df = ren_df)
    coords$label[coords$label == ""] <- coords$label_orig[coords$label == ""]
  } else {
    coords$label <- coords$label_orig
  }
  variancia <- (pcoa_result$eig / sum(pcoa_result$eig)) * 100
  list(coords = coords, var_expl = variancia, pcoa_obj = pcoa_result)
}

calcular_K_safe <- function(coords_df, tree){
  # tenta calcular K para PCoA1 e PCoA2; retorna vetor com K e p ou NA
  res <- list()
  for(i in 1:2){
    out <- tryCatch({
      trait <- coords_df[, i]
      names(trait) <- coords_df$species
      sig <- phylosig(tree, trait, method = "K", test = TRUE)
      # phylosig devolve tipicamente list com $K e $P (ou $Pval). Tentativa defensiva:
      Kval <- if(!is.null(sig$K)) sig$K else if(!is.null(sig$k)) sig$k else NA
      pval <- if(!is.null(sig$P)) sig$P else if(!is.null(sig$Pval)) sig$Pval else if(!is.null(sig$p.value)) sig$p.value else NA
      list(K = Kval, p = pval)
    }, error = function(e) list(K = NA, p = NA))
    res[[paste0("PCoA", i)]] <- out
  }
  res
}

# ---------------------------
# Plot PCoA bonito com anotações de K e p; aceita metadata (coluna 'tipo')
# ---------------------------
plot_pcoa_pretty <- function(pcoa_result, title,
                             xlim = NULL, ylim = NULL,
                             metadata = NULL, group_col = "tipo",
                             add_ellipses = TRUE,
                             point_size = 4.5,      # era 6
                             label_size = 4.5,      # era 6
                             base_size = 16,        # era 20
                             title_size = 18,       # padronizado
                             axis_title_size = 14,  # era 18
                             axis_text_size = 12){  # era 16
  coords <- pcoa_result$coords
  coords$label_expr <- coords$label  # already short code
  v1 <- round(pcoa_result$var_expl[1], 1)
  v2 <- round(pcoa_result$var_expl[2], 1)
  
  # checar metadata: se fornecida e contendo a coluna group_col, anexar
  use_meta <- FALSE
  if(!is.null(metadata) && group_col %in% colnames(metadata)){
    meta_sub <- metadata[rownames(coords), , drop = FALSE]
    if(nrow(meta_sub) == nrow(coords)) {
      coords <- cbind(coords, meta_sub)
      use_meta <- TRUE
    }
  }
  
  p <- ggplot(coords, aes(x = PCoA1, y = PCoA2)) +
    geom_hline(yintercept = 0, color = "grey90", size = 0.25) +
    geom_vline(xintercept = 0, color = "grey90", size = 0.25)
  
  if(use_meta){
    p <- p + geom_point(aes(fill = .data[[group_col]], color = .data[[group_col]]),
                        size = point_size, shape = 21, stroke = 0.6, show.legend = TRUE)
    if(add_ellipses){
      # stat_ellipse pode falhar se grupo tiver 1 obs; proteger
      p <- p + stat_ellipse(aes(color = .data[[group_col]]), level = 0.95, linetype = "dashed", show.legend = FALSE)
    }
    p <- p + scale_fill_viridis_d(option = "D") + scale_color_viridis_d(option = "D")
  } else {
    p <- p + geom_point(size = point_size, shape = 21, fill = "white", color = "black", stroke = 0.6)
  }
  
  p <- p +
    geom_text_repel(aes(label = label_expr),
                    fontface = "italic",
                    size = label_size,
                    max.overlaps = 200,
                    box.padding = 0.6) +
    labs(title = title,
         x = paste0("PCoA1 (", v1, "%)"),
         y = paste0("PCoA2 (", v2, "%)")) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    theme_minimal(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
          axis.title = element_text(face = "bold", size = axis_title_size),
          axis.text = element_text(size = axis_text_size),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey95"))
  
  # retornar gráfico e dados (para anotar K externamente)
  list(plot = p, coords = coords)
}

# ---------------------------
# Scree plot (coluna) - mas iremos colocar 3 scree em linha
# ---------------------------
scree_pcoa <- function(pcoa_result, k_max = 6,
                       title = "",
                       base_size = 14,   # era 16
                       title_size = 18,  # mesmo tamanho das PCoAs
                       text_size = 12){  # era 14
  
  use_n <- min(k_max, length(pcoa_result$var_expl))
  
  df <- data.frame(
    Axis = factor(paste0("PCoA", 1:use_n), levels = paste0("PCoA", 1:use_n)),
    Variance = pcoa_result$var_expl[1:use_n]
  )
  
  ggplot(df, aes(x = Axis, y = Variance, group = 1)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_y_continuous(
      limits  = c(0, max(df$Variance) * 1.1),
      expand  = expansion(mult = c(0.02, 0.1))
    ) +
    labs(title = title, x = NULL, y = "% explicada") +
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = title_size, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(size = text_size)
    )
}

# -------------------------
# EXECUÇÃO (garanta que dist_abaxial, dist_adaxial, dist_combined e tree estão carregados)
# -------------------------
pcoa_abaxial  <- fazer_pcoa(dist_abaxial,  k = 10, tree = tree, ren_df = renomear_Species)
pcoa_adaxial  <- fazer_pcoa(dist_adaxial,  k = 10, tree = tree, ren_df = renomear_Species)
pcoa_combined <- fazer_pcoa(dist_combined, k = 10, tree = tree, ren_df = renomear_Species)

# calcular K (defensivo)
K_abaxial  <- calcular_K_safe(pcoa_abaxial$coords, tree)
K_adaxial  <- calcular_K_safe(pcoa_adaxial$coords, tree)
K_combined <- calcular_K_safe(pcoa_combined$coords, tree)

# limites globais (mesma escala) com margem
todos_x <- c(pcoa_abaxial$coords$PCoA1, pcoa_adaxial$coords$PCoA1, pcoa_combined$coords$PCoA1)
todos_y <- c(pcoa_abaxial$coords$PCoA2, pcoa_adaxial$coords$PCoA2, pcoa_combined$coords$PCoA2)
expand_range <- function(r, mult = 0.10){ dif <- r[2] - r[1]; c(r[1] - dif * mult, r[2] + dif * mult) }
lim_x <- expand_range(range(todos_x, na.rm = TRUE))
lim_y <- expand_range(range(todos_y, na.rm = TRUE))

# detectar metadata automaticamente se presente no ambiente (nome esperado: metadata_df)
if(exists("metadata_df")){
  metadata_use <- get("metadata_df")
} else {
  metadata_use <- NULL
}

# criar PCoA plots (retornam lista com $plot)
res1 <- plot_pcoa_pretty(pcoa_abaxial,  title = "PCoA — Abaxial",  xlim = lim_x, ylim = lim_y, metadata = metadata_use, add_ellipses = TRUE)
res2 <- plot_pcoa_pretty(pcoa_adaxial,  title = "PCoA — Adaxial",  xlim = lim_x, ylim = lim_y, metadata = metadata_use, add_ellipses = TRUE)
res3 <- plot_pcoa_pretty(pcoa_combined, title = "PCoA — Combined", xlim = lim_x, ylim = lim_y, metadata = metadata_use, add_ellipses = TRUE)

p1 <- res1$plot
p2 <- res2$plot
p3 <- res3$plot

# Scree plots em linha (3)
scree_ab <- scree_pcoa(pcoa_abaxial,  k_max = 6, title = "Scree — Abaxial")
scree_ad <- scree_pcoa(pcoa_adaxial,  k_max = 6, title = "Scree — Adaxial")
scree_co <- scree_pcoa(pcoa_combined, k_max = 6, title = "Scree — Combined")

# Coluna da esquerda: PCoAs (em cima, em coluna)
left_col <- (p1 / p2 / p3) +
  plot_annotation(theme = theme(plot.title = element_text(size = 18)))

# Coluna da direita: scree (em coluna)
right_col <- (scree_ab / scree_ad / scree_co)

# Figura final: 2 colunas, 3 linhas
final_figure <- left_col | right_col

# imprimir e salvar
print(final_figure)

ggsave("Figura_PCoA_completa.png", final_figure, width = 12, height = 10, dpi = 300)
ggsave("Figura_PCoA_completa.pdf", final_figure, width = 20, height = 10)