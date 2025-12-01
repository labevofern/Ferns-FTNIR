# =============================================================================
# ANÁLISE DE COMPONENTES PRINCIPAIS
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
install.packages(c("readxl", "ggplot2", "gridExtra", "grid", "stats", "plotly"))

library(readxl)    # para importar arquivo excel
library(ggplot2)   # para plotar os gráficos
library(gridExtra) # organiza os gráficos em grade
library(grid)      # adiciona legendas e texto geral da figura em grade
library(stats)     # para fazer a PCA, porém a função prcomp que fiz a PCA ja vem na configuração base do R
library(plotly)    # para identificar os pontos no gráfico
# =============================================================================
# Iniciar análises de PCA
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar a tabela dados NIR
Matriz_NIR <- read.table(file='TABELAS/Matriz_NIR.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

#Seleciona as colunas numéricas (da oitava em diante)
dados_numericos <- Matriz_NIR[, 7:ncol(Matriz_NIR)]

#Converte as colunas para numérico
dados_numericos <- apply(dados_numericos, 2, as.numeric)

#Realiza a PCA apenas nos dados numéricos
pca_result <- prcomp(dados_numericos, scale. = TRUE)

#Realizar PCA
var_explicada <- summary(pca_result)$importance[2,]

#Para visualização
View(var_explicada)

# =============================================================================
# Esse é para fazer das espécies juntas
# =============================================================================

# Definir paleta de cores
paleta_cores <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#D37295", "#8CD17D", "#79706E"
)

# Obtém as espécies únicas na primeira coluna da matriz
especies <- unique(Matriz_NIR[, 2])

# Cria um dataframe a partir da matriz
dfespecies <- data.frame(x = pca_result$x[,1], y = pca_result$x[,2], especie = Matriz_NIR[,1])

# Multiplica os valores por 100 e formata como porcentagem
var_explicada <- var_explicada * 100

# Definir limites e quebras iguais aos do gráfico de referência
x_limits <- c(-100, 160)
y_limits <- c(-35, 35)

x_breaks <- seq(-100, 150, by = 50)
y_breaks <- seq(-35, 35, by = 10)

# Gráfico geral
g1 <- ggplot(dfespecies, aes(x = x, y = y, color = especie)) +
  geom_point(shape = 20, size = 2.5) +
  stat_ellipse(level = 0.95, type = "norm", linewidth = 0.5) +
  scale_color_manual(values = paleta_cores) +
  guides(color = guide_legend(title = NULL)) +
  labs(x = paste0("PCA1 (", sprintf("%.2f", var_explicada[1]), "%)"), 
       y = paste0("PCA2 (", sprintf("%.2f", var_explicada[2]), "%)")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 9, face = "italic"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.margin = margin(5, 5, 5, 5)
  ) +
  guides(color = guide_legend(ncol = 7, title = NULL)) +
  # USAR MESMOS LIMITES E QUEBRAS
  scale_x_continuous(limits = x_limits, breaks = x_breaks) +
  scale_y_continuous(limits = y_limits, breaks = y_breaks)

plot(g1)

# Salvar a figura combinada
ggsave("FIGURAS/Espécies_juntas.png", plot = g1, width = 8, height = 5, dpi = 300, units = "in")

# =============================================================================
# Esse é para fazer 3D e visualizar
# =============================================================================

# Cria um dataframe para o PCA com três componentes principais
dfespecies <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  Species = Matriz_NIR[, 2]
)

# Define as cores manualmente
cores_manualdfespecies <- c("DeepSkyBlue", "LawnGreen", "Gold", "MediumOrchid", "black", "Sienna", "Tomato", "DeepPink1")

# Cria o gráfico 3D
plot_ly(dfespecies, x = ~PC1, y = ~PC2, z = ~PC3, color = ~ Species, colors = cores_manualdfespecies, marker = list(size = 3)) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = paste0("PCA1 (", sprintf("%.2f", var_explicada[1] * 100), "%)")),
    yaxis = list(title = paste0("PCA2 (", sprintf("%.2f", var_explicada[2] * 100), "%)")),
    zaxis = list(title = paste0("PCA3 (", sprintf("%.2f", var_explicada[3] * 100), "%)"))
  ))

# =============================================================================
# Esse é para fazer das espécies separadas
# =============================================================================

# Gráficos individuais
plots_especies <- list()

for (i in seq_along(especies)) {
  especie_atual <- especies[i]
  cor_atual <- paleta_cores[i]
  
  df_filtrado <- dfespecies[dfespecies$especie == especie_atual, ]
  
  p <- ggplot(df_filtrado, aes(x = x, y = y)) +
    geom_point(color = cor_atual, shape = 20, size = 2.5) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),   # remove título
      axis.title = element_blank(),   # remove nomes dos eixos
      axis.text = element_text(size = 7)
    ) +
    # USAR LIMITES E QUEBRAS IDÊNTICOS AO GERAL
    scale_x_continuous(limits = x_limits, breaks = x_breaks) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks)
  
  plots_especies[[especie_atual]] <- p
}

# Combinar todos os gráficos em uma única figura
library(patchwork)

n_cols <- 5
n_rows <- ceiling(length(especies) / n_cols)

combined_plot <- wrap_plots(plots_especies, ncol = n_cols, nrow = n_rows) +
  plot_annotation(title = NULL)  # sem título geral

# Exibir
print(combined_plot)

# Salvar
ggsave("FIGURAS/Especies_individual_combinado_sem_titulos.png",
       plot = combined_plot,
       width = 12, height = 3 * n_rows, dpi = 300, units = "in")

# =============================================================================
# Abaxial e Adaxial
# =============================================================================

# Sobrescrever a coluna original, mantendo só as 3 primeiras letras
Matriz_NIR$Abaxial_or_Adaxial <- substr(Matriz_NIR$Abaxial_or_Adaxial, 1, 3)

# Agora rodar as PCAs
pca_aba <- prcomp(dados_numericos[Matriz_NIR$Abaxial_or_Adaxial == "ABA", ], scale. = TRUE)
pca_ada <- prcomp(dados_numericos[Matriz_NIR$Abaxial_or_Adaxial == "ADA", ], scale. = TRUE)

var_explicativaaba <- summary(pca_aba)$importance[2,] #Para fazer a PCA pca_aba
var_explicativaada <- summary(pca_ada)$importance[2,]

View(var_explicativaaba) #Para ver a variação especifica da PCA aba
View(var_explicativaada) #Para ver a variação especifica das PCA ada

# Determinar limites comuns para PC1 e PC2
x_lim <- range(c(pca_aba$x[, "PC1"], pca_ada$x[, "PC1"]))
y_lim <- range(c(pca_aba$x[, "PC2"], pca_ada$x[, "PC2"]))

#Multiplica os valores por 100 e formata como porcentagem
var_explicativa_aba <- var_explicativaaba * 100

# Gráfico para folhas estéreis
g5 <- ggplot(data = as.data.frame(pca_aba$x), 
             aes(x = PC1, y = PC2, 
                 color = Matriz_NIR$Specie[Matriz_NIR$`Abaxial_or_Adaxial` == "ABA"])) +
  geom_point(shape = 20, size = 2.5) +   # pontos menores e consistentes
  scale_color_manual(values = paleta_cores) +
  guides(color = guide_legend(title = NULL)) +
  labs(x = paste0("PCA1 (", sprintf("%.2f", var_explicativa_aba[1]), "%)"), 
       y = paste0("PCA2 (", sprintf("%.2f", var_explicativa_aba[2]), "%)"),
       title = "Abaxial face of fronds") +
  xlim(x_lim) + ylim(y_lim) +            # limites comuns
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),    
    axis.title = element_text(size = 15), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    legend.position = c(0.07, 0.64),  
    legend.background = element_blank(),  
    legend.text = element_text(size = 12, face = "italic", 
                               margin = margin(t = 5, r = 5, b = 5, l = 5)),  
    legend.key.height = unit(1.5, "lines"),  
    legend.key.width = unit(1.5, "lines")
  )

plot(g5)

#Multiplica os valores por 100 e formata como porcentagem
var_explicativa_ada <- var_explicativaada * 100

# Gráfico para folhas férteis
g6 <- ggplot(data = as.data.frame(pca_ada$x), 
             aes(x = PC1, y = PC2, 
                 color = Matriz_NIR$Specie[Matriz_NIR$`Abaxial_or_Adaxial` == "ADA"])) +
  geom_point(shape = 20, size = 2.5) +
  scale_color_manual(values = paleta_cores) +
  guides(color = guide_legend(title = NULL)) +
  labs(x = paste0("PCA1 (", sprintf("%.2f", var_explicativa_ada[1]), "%)"), 
       y = paste0("PCA2 (", sprintf("%.2f", var_explicativa_ada[2]), "%)"),
       title = "Adaxial face of fronds") +
  xlim(x_lim) + ylim(y_lim) + 
  theme_minimal() +
  theme(
    axis.text = element_text(size = 10),    
    axis.title = element_text(size = 15), 
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    legend.position = c(0.07, 0.64),  
    legend.background = element_blank(),  
    legend.text = element_text(size = 12, face = "italic", 
                               margin = margin(t = 5, r = 5, b = 5, l = 5)),  
    legend.key.height = unit(1.5, "lines"),  
    legend.key.width = unit(1.5, "lines")
  )

plot(g6)

ggsave("Figuras/Frondes_ABA.png", plot = g5, width = 11, height = 7, dpi= 300, units = "in")
ggsave("Figuras/Frondes_ADA.png", plot = g6, width = 11, height = 7, dpi= 300, units = "in")
