# =============================================================================
# VER ESPECTROS E OUTLINES NIR
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
install.packages(c("readr", "readxl", "ggplot2", "dplyr"))

library(readr)   # Leitura eficiente de dados (CSV, TXT)
library(readxl)  # Leitura de arquivos Excel (XLSX, XLS)
library(ggplot2) # Visualização de dados avançada
library(dplyr)   # Manipulação de dados (filtros, seleção, agrupamento)
# =============================================================================
# CHAMAR A TABELA
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar a tabela dados NIR
Matriz_NIR <- read.table(file='TABELAS/Matriz_NIR.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

# =============================================================================
# Opção 1 - espectros médios das espécies
# =============================================================================

# Identificar colunas espectrais
nir_cols <- grep("^X", names(Matriz_NIR), value = TRUE)

# Converter colunas NIR para numérico
Matriz_NIR[nir_cols] <- lapply(Matriz_NIR[nir_cols], function(x) as.numeric(as.character(x)))

# Preparar dados para o gráfico (long format)
dados_grafico <- Matriz_NIR %>%
  tidyr::pivot_longer(
    cols = all_of(nir_cols), 
    names_to = "Wavelength", 
    values_to = "Absorbance"
  ) %>%
  mutate(Wavelength = as.numeric(sub("^X", "", Wavelength)))

# Calcular a média de absorbância por espécie e comprimento de onda
dados_media <- dados_grafico %>%
  group_by(Specie, Wavelength) %>%
  summarise(Absorbance = mean(Absorbance, na.rm = TRUE), .groups = "drop")

# Definir limites manuais do eixo Y
y_lim_inf <- 0.2
y_lim_sup <- 0.8

paleta_cores <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
  "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#D37295", "#8CD17D", "#79706E"
)

# Criar gráfico com médias - COM LIMITES Y DEFINIDOS
p1 <- ggplot(dados_media, aes(x = Wavelength, y = Absorbance, color = Specie)) +
  geom_line(size = 0.5) +
  scale_color_manual(values = paleta_cores, name = "Species") +
  scale_y_continuous(limits = c(y_lim_inf, y_lim_sup),
                     breaks = c(0.2, 0.4, 0.6, 0.8)) +  # Limites e quebras definidos
  labs(
    title = "",
    x = expression("Wavenumber (cm"^{-1}*")"),
    y = "Absorbance"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.82, 0.70),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  guides(color = guide_legend(ncol = 2))

print(p1)

ggsave("FIGURAS/Espectros_médios_complicadas.png", plot = p1, width = 8, height = 5, dpi = 300)

# =============================================================================
# Opção 1 - espectros médios das espécies (com desvio padrão)
# =============================================================================

# Calcular a média e desvio padrão por espécie e comprimento de onda
dados_estatisticas <- dados_grafico %>%
  group_by(Specie, Wavelength) %>%
  summarise(
    Media = mean(Absorbance, na.rm = TRUE),
    Desvio_Padrao = sd(Absorbance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Limite_Superior = Media + Desvio_Padrao,
    Limite_Inferior = Media - Desvio_Padrao
  )

# Criar gráfico com médias e desvio padrão - MESMOS LIMITES Y
dv <- ggplot(dados_estatisticas, aes(x = Wavelength, y = Media, color = Specie)) +
  geom_line(size = 0.5) +
  geom_ribbon(aes(ymin = Limite_Inferior, ymax = Limite_Superior, fill = Specie),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = paleta_cores, name = "Species") +
  scale_fill_manual(values = paleta_cores, name = "Species") +
  scale_y_continuous(limits = c(y_lim_inf, y_lim_sup),
                     breaks = c(0.2, 0.4, 0.6, 0.8)) +  # MESMOS LIMITES E QUEBRAS
  labs(
    title = "",
    x = expression("Wavenumber (cm"^{-1}*")"),
    y = "Absorbance"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.83, 0.70),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  guides(color = guide_legend(ncol = 2))

print(dv)

ggsave("FIGURAS/Espectros_médios_DESVIO_PADRÃO_complicadas.png", plot = dv, width = 8, height = 5, dpi = 300)

# =============================================================================
# Opção 2 - Todos os espectros de todos os espécimes
# =============================================================================

# Identificar colunas espectrais
nir_cols <- grep("^X", names(Matriz_NIR), value = TRUE)

# Converter colunas NIR para numérico
Matriz_NIR[nir_cols] <- lapply(Matriz_NIR[nir_cols], function(x) as.numeric(as.character(x)))

# Adicionar um ID único para cada linha (cada espectro individual)
Matriz_NIR <- Matriz_NIR %>%
  mutate(Espectro_ID = row_number())

# Preparar dados para o gráfico (long format) - APENAS DADOS BRUTOS
dados_grafico <- Matriz_NIR %>%
  tidyr::pivot_longer(
    cols = all_of(nir_cols), 
    names_to = "Wavelength", 
    values_to = "Absorbance"
  ) %>%
  mutate(Wavelength = as.numeric(sub("^X", "", Wavelength)))

# Criar gráfico APENAS com espectros brutos individuais (SEM LIMITES Y)
p1 <- ggplot(dados_grafico, aes(x = Wavelength, y = Absorbance, color = Specie)) +
  geom_line(aes(group = Espectro_ID),  # Usando o ID único criado
            size = 0.3, alpha = 0.7) +
  scale_color_manual(values = paleta_cores, name = "Species") +
  labs(
    title = "",
    x = expression("Wavenumber (cm"^{-1}*")"),
    y = "Absorbance"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.82, 0.70),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  ) +
  guides(color = guide_legend(ncol = 2))

print(p1)

ggsave("FIGURAS/Espectros_bruto_complicadas.png", plot = p1, width = 8, height = 5, dpi = 300)

# =============================================================================
# VISUALIZAR ESPECTRO POR AMOSTRA
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar a tabela média dos dados NIR
Matriz_NIR <- read.table(file='TABELAS/Matriz_NIR.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

categoria = "Collector_and_number"  # nome da coluna que identifica indivíduos
unique(Matriz_NIR$Collector_and_number)

AbaAda = 'Abaxial_or_Adaxial'  # nome da coluna com categoria para colorir espectros

# Extrair apenas as 3 primeiras letras da coluna Abaxial_or_Adaxial
Matriz_NIR$Abaxial_or_Adaxial <- substr(Matriz_NIR$Abaxial_or_Adaxial, 1, 3)

# Verificar os valores únicos após a modificação
unique(Matriz_NIR$Abaxial_or_Adaxial)

# Pega colunas NIR (opcional, apenas se for usar depois)
cls = grep("X", colnames(Matriz_NIR), ignore.case = F)

# Individuos unicos
indv = Matriz_NIR[, categoria]
uindv = unique(indv)

cat = Matriz_NIR[, AbaAda]                  
cores = c("red","blue") # Para abaxial e adaxial
cat.cor = cores[as.numeric(as.factor(cat))]
names(cat.cor) = cat

# Gera um pdf
pdf(file = 'FIGURAS/Espectros por espécime.pdf', width = 8.5, height = 10)
 
# Divide a página em tres linhas e duas colunas (seis graficos por pagina)
par(mfrow=c(2,2))

# Para cada individuo unico plota os spectros colorindo abaxial e adaxial
u = 1
for (u in 1:length(uindv)) {
  # Filtra os dados para o indivíduo
  d = Matriz_NIR[indv == uindv[u], cls]
  
  faccore = cat.cor[indv == uindv[u]]
  
  # Pega eixo x do nome das variáveis NIR	
  xx = colnames(d)
  xx = gsub("X", "", xx)
  xx = as.numeric(xx)
  
  # Define a faixa desejada de absorbância
  yl = c(0.2, 0.9)
  
  # Plot uma figura vazia
  plot(xx, d[1,], type = 'n', xlab = 'expression("Wavenumber (cm"^{-1}*")")', ylab = 'Absorbance', ylim = yl)
  
  # Plot cada espectro
  for (n in 1:nrow(d)) {
    y = d[n,]
    # Ajusta os valores de absorbância para a faixa desejada
    y_adj = pmin(pmax(y, 0.2), 0.8)
    points(xx, y_adj, type = 'l', col = faccore[n])
  }
  
  legend("topleft", legend = uindv[u], bty = 'n', inset = 0.0)
  legend("topright", legend = unique(names(cat.cor)), bty = 'n', inset = 0.0, lwd = 2, col = unique(cat.cor))
}

dev.off()

# =============================================================================
# VISUALIZAR ESPECTRO POR ESPÉCIE
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar a tabela média dos dados NIR
Matriz_NIR <- read.table(file='TABELAS/Matriz_NIR.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

categoria = "Specie"  # nome da coluna que identifica indivíduos
unique(Matriz_NIR$Specie)

AbaAda = 'Abaxial_or_Adaxial'  # nome da coluna com categoria para colorir espectros

# Extrair apenas as 3 primeiras letras da coluna Abaxial_or_Adaxial
Matriz_NIR$Abaxial_or_Adaxial <- substr(Matriz_NIR$Abaxial_or_Adaxial, 1, 3)

# Verificar os valores únicos após a modificação
unique(Matriz_NIR$Abaxial_or_Adaxial)

# Pega colunas NIR (opcional, apenas se for usar depois)
cls = grep("X", colnames(Matriz_NIR), ignore.case = F)

# Individuos unicos
indv = Matriz_NIR[, categoria]
uindv = unique(indv)

cat = Matriz_NIR[, AbaAda]                  
cores = c("red","blue") # Para abaxial e adaxial
cat.cor = cores[as.numeric(as.factor(cat))]
names(cat.cor) = cat

# Gera um pdf
pdf(file = 'FIGURAS/teste.pdf', width = 8.5, height = 10)

# Divide a página em tres linhas e duas colunas (seis graficos por pagina)
par(mfrow=c(2,2))

# Para cada individuo unico plota os spectros colorindo abaxial e adaxial
u = 1
for (u in 1:length(uindv)) {
  # Filtra os dados para o indivíduo
  d = Matriz_NIR[indv == uindv[u], cls]
  
  faccore = cat.cor[indv == uindv[u]]
  
  # Pega eixo x do nome das variáveis NIR	
  xx = colnames(d)
  xx = gsub("X", "", xx)
  xx = as.numeric(xx)
  
  # Define a faixa desejada de absorbância
  yl = c(0.2, 0.9)
  
  # Plot uma figura vazia
  plot(xx, d[1,], type = 'n', xlab = 'expression("Wavenumber (cm"^{-1}*")")', ylab = 'Absorbance', ylim = yl)
  
  # Plot cada espectro
  for (n in 1:nrow(d)) {
    y = d[n,]
    # Ajusta os valores de absorbância para a faixa desejada
    y_adj = pmin(pmax(y, 0.2), 0.8)
    points(xx, y_adj, type = 'l', col = faccore[n])
  }
  
  legend("topleft", legend = uindv[u], bty = 'n', inset = 0.0)
  legend("topright", legend = unique(names(cat.cor)), bty = 'n', inset = 0.0, lwd = 2, col = unique(cat.cor))
}
dev.off()