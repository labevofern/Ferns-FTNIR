# =============================================================================
# PERMANOVA E TESTE DE DISPERSÃO
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
install.packages(c("vegan", "ggplot2", "dplyr", "pheatmap", "RColorBrewer", "viridis", "combinat"))

library(vegan)          # Análise ecológica e de comunidade: diversidade, ordenação
library(ggplot2)        # Criação de gráficos elegantes e personalizáveis
library(dplyr)          # Manipulação de dados: filtrar, selecionar, modificar
library(pheatmap)       # Criação de heatmaps/tabelas de calor com anotações
library(RColorBrewer)   # Paletas de cores atrativas para visualizações
library(viridis)        # Paletas de cores colorblind-friendly e perceptualmente uniformes
library(combinat)       # Combinatória: permutações, combinações e outras operações
library(reshape2)       # 
# =============================================================================
# Importar tabela
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar a tabela dados NIR
dados <- read.table(file='TABELAS/Matriz_NIR.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

# Primeiro, converter colunas espectrais para numérico (colunas 7 até o final)
dados[, 7:ncol(dados)] <- lapply(dados[, 7:ncol(dados)], function(x) as.numeric(as.character(x)))

# ==============================================================================
# Calcular a média por indivíduo
# ==============================================================================

# Calcular médias por Collector_and_number para as colunas espectrais
medias <- dados %>%
  group_by(Collector_and_number) %>%
  summarise(across(6:last_col(), mean, na.rm = TRUE))

# Obter as informações únicas de Specie para cada Collector_and_number
especies_unicas <- dados %>%
  distinct(Collector_and_number, Specie)

# Combinar os resultados
dados <- medias %>%
  left_join(especies_unicas, by = "Collector_and_number") %>%
  select(Specie, Collector_and_number, everything())

# Visualizar o resultado
print(dados)

# ==============================================================================
# Permanova e dispersão
# ==============================================================================

dados_espectros <- dados[, 3:ncol(dados)]  # colunas 7 em diante
species <- dados$Specie                   # vetor de espécies

dados_espectros <- dados_espectros %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

# Transformação log
# Apenas escala os dados sem log
dados_scaled <- scale(dados_espectros)

# Matriz de distância
dist_matrix <- vegdist(dados_scaled, method = "euclidean")

# PERMANOVA
adonis_result <- adonis2(dist_matrix ~ species, data = dados, permutations = 9999)
print(adonis_result)

# Teste de dispersão
disp <- betadisper(dist_matrix, species)
anova_disp <- anova(disp)
print(anova_disp)

# Transformar objeto betadisper em data frame
disp_df <- data.frame(Species = species, Distance = disp$distances)

# Definir paleta de cores
paleta_cores <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#D37295", "#8CD17D", "#79706E"
)

# Calcular os centroides por espécie
centroides <- disp_df %>%
  group_by(Species) %>%
  summarise(Centroide = mean(Distance))

g1 <- ggplot(disp_df, aes(x = Species, y = Distance, fill = Species)) +
  geom_boxplot() +
  geom_point(data = centroides, aes(x = Species, y = Centroide),
             color = "black", size = 3, shape = 20) +
  scale_fill_manual(values = paleta_cores) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none"  # remove a legenda
  ) +
  labs(title = "Intraspecies spectral variability in relation to the centroid",
       y = "Distance to centroid", x = "")

plot(g1)

# Salvar a figura combinada
ggsave("FIGURAS/Variabilidade espectral.png", plot = g1, width = 8, height = 5, dpi = 300, units = "in")

# ==============================================================================
# Permanova par a par
# ==============================================================================

# Obter nomes únicos das espécies
esp <- unique(species)

# Criar todas as combinações de pares de espécies
pares <- combn(esp, 2, simplify = FALSE)

results_permanova <- data.frame(Especies = character(),
                                F = numeric(),
                                R2 = numeric(),
                                p = numeric(),
                                stringsAsFactors = FALSE)

for (p in pares) {
  # Subset dos dados só para o par
  idx <- species %in% p
  dist_sub <- vegdist(dados_scaled[idx, ], method = "euclidean")
  species_sub <- species[idx]
  
  # PERMANOVA
  adonis_sub <- adonis2(dist_sub ~ species_sub, permutations = 999)
  
  results_permanova <- rbind(results_permanova,
                             data.frame(Especies = paste(p, collapse = " vs "),
                                        F = adonis_sub$F[1],
                                        R2 = adonis_sub$R2[1],
                                        p = adonis_sub$`Pr(>F)`[1]))
}

print(results_permanova)

# ==============================================================================
# Dispersão par a par
# ==============================================================================

results_disp <- data.frame(Especies = character(),
                           F = numeric(),
                           p = numeric(),
                           stringsAsFactors = FALSE)

for (p in pares) {
  # Subset dos dados só para o par
  idx <- species %in% p
  dist_sub <- vegdist(dados_scaled[idx, ], method = "euclidean")
  species_sub <- species[idx]
  
  # Betadisper (dispersão)
  disp_sub <- betadisper(dist_sub, species_sub)
  anova_sub <- anova(disp_sub)
  
  results_disp <- rbind(results_disp,
                        data.frame(Especies = paste(p, collapse = " vs "),
                                   F = anova_sub$`F value`[1],
                                   p = anova_sub$`Pr(>F)`[1]))
}

print(results_disp)

# ==============================================================================
# Gráfico de heatmap
# ==============================================================================
# Carregar bibliotecas necessárias
library(ggplot2)
library(reshape2)
library(dplyr)

# Criar dataframe com os dados
data <- data.frame(
  Comparison = c(
    "ADIPUL vs MENSER", "ADIPUL vs SERTRI", "ADIPUL vs CYCMEN", "ADIPUL vs BLEOCC",
    "ADIPUL vs LYGVEN", "ADIPUL vs DIPCRI", "ADIPUL vs CYAMIC", "ADIPUL vs MICGEM",
    "ADIPUL vs NIPCRA", "ADIPUL vs PLEAST", "ADIPUL vs ANEHIR", "ADIPUL vs DICFLE",
    "MENSER vs SERTRI", "MENSER vs CYCMEN", "MENSER vs BLEOCC", "MENSER vs LYGVEN",
    "MENSER vs DIPCRI", "MENSER vs CYAMIC", "MENSER vs MICGEM", "MENSER vs NIPCRA",
    "MENSER vs PLEAST", "MENSER vs ANEHIR", "MENSER vs DICFLE", "SERTRI vs CYCMEN",
    "SERTRI vs BLEOCC", "SERTRI vs LYGVEN", "SERTRI vs DIPCRI", "SERTRI vs CYAMIC",
    "SERTRI vs MICGEM", "SERTRI vs NIPCRA", "SERTRI vs PLEAST", "SERTRI vs ANEHIR",
    "SERTRI vs DICFLE", "CYCMEN vs BLEOCC", "CYCMEN vs LYGVEN", "CYCMEN vs DIPCRI",
    "CYCMEN vs CYAMIC", "CYCMEN vs MICGEM", "CYCMEN vs NIPCRA", "CYCMEN vs PLEAST",
    "CYCMEN vs ANEHIR", "CYCMEN vs DICFLE", "BLEOCC vs LYGVEN", "BLEOCC vs DIPCRI",
    "BLEOCC vs CYAMIC", "BLEOCC vs MICGEM", "BLEOCC vs NIPCRA", "BLEOCC vs PLEAST",
    "BLEOCC vs ANEHIR", "BLEOCC vs DICFLE", "LYGVEN vs DIPCRI", "LYGVEN vs CYAMIC",
    "LYGVEN vs MICGEM", "LYGVEN vs NIPCRA", "LYGVEN vs PLEAST", "LYGVEN vs ANEHIR",
    "LYGVEN vs DICFLE", "DIPCRI vs CYAMIC", "DIPCRI vs MICGEM", "DIPCRI vs NIPCRA",
    "DIPCRI vs PLEAST", "DIPCRI vs ANEHIR", "DIPCRI vs DICFLE", "CYAMIC vs MICGEM",
    "CYAMIC vs NIPCRA", "CYAMIC vs PLEAST", "CYAMIC vs ANEHIR", "CYAMIC vs DICFLE",
    "MICGEM vs NIPCRA", "MICGEM vs PLEAST", "MICGEM vs ANEHIR", "MICGEM vs DICFLE",
    "NIPCRA vs PLEAST", "NIPCRA vs ANEHIR", "NIPCRA vs DICFLE", "PLEAST vs ANEHIR",
    "PLEAST vs DICFLE", "ANEHIR vs DICFLE"
  ),
  F = c(
    17.77564, 7.371135, 51.02177, 11.82235, 3.200064, 36.19, 2.284757, 5.37933,
    9.052041, 9.979181, 21.28268, 6.052548, 0.222507, 4.336854, 1.417138, 3.802675,
    0.481709, 23.46683, 12.23269, 19.99019, 27.18183, 1.116202, 32.9251, 2.140656,
    0.393983, 1.645237, 0.773056, 10.25755, 5.657112, 12.58722, 16.65455, 1.219777,
    13.57794, 11.24512, 13.94203, 8.795429, 48.63422, 18.07814, 20.3703, 30.55774,
    9.757097, 75.19284, 1.564991, 3.842676, 16.61741, 8.368707, 16.02651, 22.07394,
    3.000392, 26.70743, 7.618167, 6.956593, 5.592464, 12.31401, 15.18267, 4.518898,
    11.78411, 42.23029, 19.03784, 24.89468, 35.07648, 1.640265, 66.12787, 2.554868,
    5.263056, 5.784821, 29.71349, 1.516314, 3.512368, 5.945015, 17.02959, 4.254509,
    0.855412, 23.51166, 3.059748, 30.8544, 2.335609, 41.67172
  ),
  R2 = c(
    0.496864, 0.302453, 0.772802, 0.396426, 0.150946, 0.667835, 0.107342, 0.230089,
    0.334616, 0.356665, 0.541783, 0.27446, 0.01292, 0.224279, 0.072984, 0.174413,
    0.026064, 0.552592, 0.404618, 0.526193, 0.60161, 0.05839, 0.672969, 0.132625,
    0.022651, 0.088239, 0.043496, 0.363002, 0.249684, 0.425428, 0.494868, 0.066948,
    0.47512, 0.428465, 0.481723, 0.369627, 0.752453, 0.546528, 0.575915, 0.670747,
    0.394113, 0.852596, 0.079989, 0.175925, 0.466553, 0.317373, 0.471001, 0.55083,
    0.142873, 0.625358, 0.297374, 0.268009, 0.237045, 0.406215, 0.457548, 0.200671,
    0.424131, 0.689696, 0.514011, 0.580368, 0.660867, 0.083515, 0.805182, 0.118529,
    0.216916, 0.233402, 0.609964, 0.081891, 0.163272, 0.248278, 0.486149, 0.210052,
    0.045367, 0.566387, 0.160535, 0.631558, 0.127381, 0.722568
  ),
  p = c(
    0.001, 0.01, 0.001, 0.002, 0.089, 0.001, 0.13, 0.011, 0.003, 0.004, 0.002, 0.017,
    0.691, 0.026, 0.25, 0.057, 0.581, 0.001, 0.003, 0.001, 0.001, 0.305, 0.002, 0.212,
    0.567, 0.21, 0.408, 0.004, 0.023, 0.005, 0.001, 0.275, 0.003, 0.002, 0.003, 0.003,
    0.001, 0.001, 0.002, 0.001, 0.001, 0.001, 0.216, 0.043, 0.001, 0.007, 0.001, 0.001,
    0.078, 0.001, 0.007, 0.012, 0.029, 0.001, 0.001, 0.05, 0.007, 0.001, 0.001, 0.001,
    0.001, 0.197, 0.001, 0.098, 0.029, 0.022, 0.001, 0.223, 0.064, 0.022, 0.001, 0.038,
    0.363, 0.001, 0.089, 0.001, 0.129, 0.002
  )
)

# Extrair todas as espécies únicas
all_species <- unique(c(
  sapply(strsplit(data$Comparison, " vs "), function(x) x[1]),
  sapply(strsplit(data$Comparison, " vs "), function(x) x[2])
))

# Criar matriz vazia
n_species <- length(all_species)
p_matrix <- matrix(NA, nrow = n_species, ncol = n_species,
                   dimnames = list(all_species, all_species))

# Preencher a diagonal com 1 (auto-comparações)
diag(p_matrix) <- 1

# Preencher a matriz com valores-p
for (i in 1:nrow(data)) {
  comp <- strsplit(data$Comparison[i], " vs ")[[1]]
  sp1 <- comp[1]
  sp2 <- comp[2]
  p_val <- data$p[i]
  
  p_matrix[sp1, sp2] <- p_val
  p_matrix[sp2, sp1] <- p_val
}

# Converter para formato longo para ggplot
p_df <- melt(p_matrix)
colnames(p_df) <- c("Species1", "Species2", "p_value")

# Plotar Dotplot
g1 = ggplot(p_df, aes(x = Species1, y = Species2)) +
  geom_point(aes(size = -log10(p_value), color = p_value < 0.05), alpha = 0.8) +
  scale_size_continuous(name = "-log10(p)", range = c(1, 8)) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey70"),
                     name = "Significant") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),  # Remove título do eixo x
    axis.title.y = element_blank()   # Remove título do eixo y
  ) +
  coord_fixed() +
  ggtitle("") +
  labs(x = NULL, y = NULL)  # Também remove os títulos dos eixos

# Exibir o gráfico
print(g1)

ggsave("FIGURAS/Dotplot com valores de p.png", width = 7, height = 7, dpi = 300)
