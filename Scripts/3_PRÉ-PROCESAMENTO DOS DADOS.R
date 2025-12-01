# =============================================================================
# PRÉ-PROCESAMENTO DE DADOS
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
install.packages(c("ggplot2", "tidyr", "dplyr", "plorly"))

library(ggplot2)  # Visualização de dados avançada
library(tidyr)    # Manipulação de dados (pivotagem, nesting, etc.)
library(dplyr)    # Manipulação de dados (filtros, seleção, agrupamento)
library(plotly)   # Para identificar os pontos no gráfico
# =============================================================================
# Standard Normal Variate (SNV)
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar as tabelas NIR com os dados médios
Matriz_NIR <- read.table(file='TABELAS/Matriz_NIR.csv', header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

# Função para aplicar SNV a um conjunto de dados
processar_dados <- function(dados, nome_conjunto) {
  # Função para aplicar SNV a um único espectro
  apply_snv <- function(spectrum) {
    mean_spectrum <- mean(as.numeric(spectrum))
    sd_spectrum <- sd(as.numeric(spectrum))
    snv_spectrum <- (as.numeric(spectrum) - mean_spectrum) / sd_spectrum
    return(snv_spectrum)
  }
  
  # Separar os dados em informações categóricas e valores espectrais
  dados_categoricos <- dados[, 1:6]
  valores_espectrais <- dados[, 7:ncol(dados)]
  
  # Garantir que valores_espectrais contenha apenas números
  valores_espectrais <- apply(valores_espectrais, 2, as.numeric)
  
  # Aplicar SNV aos valores espectrais
  valores_espectrais_snv <- t(apply(valores_espectrais, 1, apply_snv))
  
  # Manter os nomes originais das colunas
  colnames(valores_espectrais_snv) <- colnames(valores_espectrais)
  
  # Combinar de volta as informações categóricas com os valores espectrais normalizados
  data_snv <- cbind(dados_categoricos, valores_espectrais_snv)
  
  # Identificar colunas espectrais
  nir_cols <- grep("^X", names(data_snv), value = TRUE)
  
  # Converter colunas NIR para numérico
  data_snv[nir_cols] <- lapply(data_snv[nir_cols], function(x) as.numeric(as.character(x)))
  
  # Adicionar um ID único para cada linha (cada espectro individual)
  data_snv <- data_snv %>%
    mutate(Espectro_ID = row_number())
  
  # Preparar dados para o gráfico
  dados_grafico <- data_snv %>%
    tidyr::pivot_longer(cols = all_of(nir_cols), 
                        names_to = "Wavelength", 
                        values_to = "Absorbance") %>%
    mutate(Wavelength = as.numeric(sub("^X", "", Wavelength)))
  
  # Definir paleta de cores
  paleta_cores <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#D37295", "#8CD17D", "#79706E"
  )
  
  # Criar o gráfico com as configurações solicitadas
  p <- ggplot(dados_grafico, aes(x = Wavelength, y = Absorbance, color = Specie)) +
    geom_line(aes(group = Espectro_ID), size = 0.3, alpha = 0.7) +
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
  
  # Salvar arquivo
  write.table(data_snv,
              file = paste0("TABELAS/Matriz_NIR_snv", nome_conjunto, ".csv"),
              row.names = FALSE,
              sep = ";",
              quote = FALSE,
              dec = ".",
              fileEncoding = "UTF-8")
  
  # Salvar figura
  ggsave(paste0("FIGURAS/Espectros processados_snv", nome_conjunto, ".png"),
         plot = p, width = 8, height = 5, dpi = 300, units = "in")
  
  # Gráfico interativo com plotly
  p_interativo <- ggplot(dados_grafico, aes(x = Wavelength, y = Absorbance,
                                            color = Specie,
                                            group = Espectro_ID,
                                            text = paste("Specie:", Specie,
                                                         "<br>Collector:", Collector_and_number,
                                                         "<br>Surface:", Abaxial_or_Adaxial))) +
    geom_line(alpha = 0.7, size = 0.3) +
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
  
  # Converter para plotly
  plotly_plot <- ggplotly(p_interativo, tooltip = "text")
  
  # Retornar os dados processados e o gráfico (opcional)
  return(list(data = data_snv, plot = p, interactive_plot = plotly_plot))
}

# Processar cada conjunto de dados
resultados <- processar_dados(Matriz_NIR, "NIR")

# Mostrar os gráficos (opcional)
print(resultados$plot)

# Mostrar gráficos interativos
resultados$interactive_plot

htmlwidgets::saveWidget(resultados$interactive_plot, "FIGURAS/Spectra_interactive_AB_3001.html")












# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Carregar pacotes necessários
library(ggplot2)
library(plotly)
library(tidyr)
library(dplyr)
library(prospectr) # Para derivadas e Savitzky-Golay

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar as tabelas NIR com os dados médios
Matriz_NIR <- read.table(file='TABELAS/Matriz_NIR.csv', header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

# Função para aplicar diferentes pré-processamentos
aplicar_preprocessamento <- function(dados, metodo) {
  
  # Separar os dados em informações categóricas e valores espectrais
  dados_categoricos <- dados[, 1:6]
  valores_espectrais <- dados[, 7:ncol(dados)]
  
  # Garantir que valores_espectrais contenha apenas números
  valores_espectrais <- apply(valores_espectrais, 2, as.numeric)
  comprimento_original <- ncol(valores_espectrais)
  nomes_originais <- colnames(valores_espectrais)
  
  # Aplicar diferentes métodos de pré-processamento
  if (metodo == "bruto") {
    # Dados brutos - apenas converter para numérico
    valores_processados <- valores_espectrais
    
  } else if (metodo == "snv") {
    # Aplicar SNV
    apply_snv <- function(spectrum) {
      mean_spectrum <- mean(as.numeric(spectrum))
      sd_spectrum <- sd(as.numeric(spectrum))
      snv_spectrum <- (as.numeric(spectrum) - mean_spectrum) / sd_spectrum
      return(snv_spectrum)
    }
    valores_processados <- t(apply(valores_espectrais, 1, apply_snv))
    
  } else if (metodo == "derivada1") {
    # Primeira derivada (Savitzky-Golay)
    valores_processados <- savitzkyGolay(X = valores_espectrais, m = 1, p = 2, w = 11)
    # Ajustar nomes das colunas para manter o mesmo comprimento
    colnames(valores_processados) <- nomes_originais[6:(comprimento_original-5)]
    
  } else if (metodo == "derivada2") {
    # Segunda derivada (Savitzky-Golay)
    valores_processados <- savitzkyGolay(X = valores_espectrais, m = 2, p = 2, w = 11)
    colnames(valores_processados) <- nomes_originais[6:(comprimento_original-5)]
    
  } else if (metodo == "snv_derivada1") {
    # SNV seguido de primeira derivada
    apply_snv <- function(spectrum) {
      mean_spectrum <- mean(as.numeric(spectrum))
      sd_spectrum <- sd(as.numeric(spectrum))
      snv_spectrum <- (as.numeric(spectrum) - mean_spectrum) / sd_spectrum
      return(snv_spectrum)
    }
    snv_data <- t(apply(valores_espectrais, 1, apply_snv))
    valores_processados <- savitzkyGolay(X = snv_data, m = 1, p = 2, w = 11)
    colnames(valores_processados) <- nomes_originais[6:(comprimento_original-5)]
    
  } else if (metodo == "derivada1_snv") {
    # Primeira derivada seguida de SNV
    derivada_data <- savitzkyGolay(X = valores_espectrais, m = 1, p = 2, w = 11)
    colnames(derivada_data) <- nomes_originais[6:(comprimento_original-5)]
    
    apply_snv <- function(spectrum) {
      mean_spectrum <- mean(as.numeric(spectrum))
      sd_spectrum <- sd(as.numeric(spectrum))
      snv_spectrum <- (as.numeric(spectrum) - mean_spectrum) / sd_spectrum
      return(snv_spectrum)
    }
    valores_processados <- t(apply(derivada_data, 1, apply_snv))
    
  } else if (metodo == "snv_derivada2") {
    # SNV seguido de segunda derivada
    apply_snv <- function(spectrum) {
      mean_spectrum <- mean(as.numeric(spectrum))
      sd_spectrum <- sd(as.numeric(spectrum))
      snv_spectrum <- (as.numeric(spectrum) - mean_spectrum) / sd_spectrum
      return(snv_spectrum)
    }
    snv_data <- t(apply(valores_espectrais, 1, apply_snv))
    valores_processados <- savitzkyGolay(X = snv_data, m = 2, p = 2, w = 11)
    colnames(valores_processados) <- nomes_originais[6:(comprimento_original-5)]
    
  } else if (metodo == "golay") {
    # Filtro Savitzky-Golay (suavização)
    valores_processados <- savitzkyGolay(X = valores_espectrais, m = 0, p = 2, w = 11)
    colnames(valores_processados) <- nomes_originais[6:(comprimento_original-5)]
  }
  
  # Combinar de volta as informações categóricas com os valores processados
  if (metodo %in% c("derivada1", "derivada2", "snv_derivada1", "derivada1_snv", "snv_derivada2", "golay")) {
    # Para métodos que reduzem o número de colunas, ajustar os dados categóricos
    dados_processados <- cbind(dados_categoricos, valores_processados)
  } else {
    dados_processados <- cbind(dados_categoricos, valores_processados)
  }
  
  return(dados_processados)
}

# Função para criar e salvar gráficos
criar_grafico_espectros <- function(dados_processados, nome_metodo) {
  
  # Identificar colunas espectrais
  nir_cols <- grep("^X", names(dados_processados), value = TRUE)
  
  # Converter colunas NIR para numérico
  dados_processados[nir_cols] <- lapply(dados_processados[nir_cols], function(x) as.numeric(as.character(x)))
  
  # Adicionar um ID único para cada linha
  dados_processados <- dados_processados %>%
    mutate(Espectro_ID = row_number())
  
  # Preparar dados para o gráfico
  dados_grafico <- dados_processados %>%
    tidyr::pivot_longer(cols = all_of(nir_cols), 
                        names_to = "Wavelength", 
                        values_to = "Absorbance") %>%
    mutate(Wavelength = as.numeric(sub("^X", "", Wavelength)))
  
  # Definir paleta de cores
  paleta_cores <- c(
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#D37295", "#8CD17D", "#79706E"
  )
  
  # Criar o gráfico
  p <- ggplot(dados_grafico, aes(x = Wavelength, y = Absorbance, color = Specie)) +
    geom_line(aes(group = Espectro_ID), linewidth = 0.3, alpha = 0.7) +
    scale_color_manual(values = paleta_cores, name = "Species") +
    labs(
      title = paste("Método:", nome_metodo),
      x = expression("Wavenumber (cm"^{-1}*")"),
      y = "Absorbance"
    ) +
    theme_minimal() +
    theme(
      legend.position.inside = c(0.82, 0.70),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    guides(color = guide_legend(ncol = 2))
  
  # Salvar figura
  ggsave(paste0("FIGURAS/Espectros_", nome_metodo, ".png"),
         plot = p, width = 8, height = 5, dpi = 300, units = "in")
  
  # Salvar dados processados
  write.table(dados_processados,
              file = paste0("TABELAS/Matriz_NIR_", nome_metodo, ".csv"),
              row.names = FALSE,
              sep = ";",
              quote = FALSE,
              dec = ".",
              fileEncoding = "UTF-8")
  
  return(p)
}

# Lista de métodos a testar (vamos começar com os mais simples)
metodos <- c("bruto", "snv", "derivada1", "golay")

nomes_metodos <- c("Dados_Brutos", "SNV", "Primeira_Derivada", "Savitzky_Golay")

# Aplicar todos os métodos e criar gráficos
graficos <- list()

for (i in seq_along(metodos)) {
  cat("Processando:", nomes_metodos[i], "\n")
  
  tryCatch({
    # Aplicar pré-processamento
    dados_processados <- aplicar_preprocessamento(Matriz_NIR, metodos[i])
    
    # Criar e salvar gráfico
    grafico <- criar_grafico_espectros(dados_processados, nomes_metodos[i])
    
    # Armazenar gráfico na lista
    graficos[[nomes_metodos[i]]] <- grafico
    
    cat("Sucesso:", nomes_metodos[i], "\n")
    
  }, error = function(e) {
    cat("Erro no método", nomes_metodos[i], ":", e$message, "\n")
  })
}

# Mostrar todos os gráficos (opcional)
for (nome in names(graficos)) {
  print(graficos[[nome]])
}

# Criar um gráfico comparativo simples
cat("Criando gráfico comparativo...\n")

# Combinar dados dos métodos que funcionaram
dados_comparativos <- data.frame()
metodos_funcionando <- names(graficos)

for (metodo in metodos_funcionando) {
  try({
    dados_metodo <- read.table(paste0("TABELAS/Matriz_NIR_", metodo, ".csv"),
                               header = TRUE, sep = ";", dec = ".")
    
    # Selecionar algumas amostras representativas
    amostras_representativas <- dados_metodo %>%
      group_by(Specie) %>%
      slice(1) %>%  # Pegar 1 amostra de cada espécie
      ungroup() %>%
      mutate(Metodo = metodo)
    
    dados_comparativos <- bind_rows(dados_comparativos, amostras_representativas)
  })
}

# Preparar dados para gráfico comparativo
if (nrow(dados_comparativos) > 0) {
  nir_cols <- grep("^X", names(dados_comparativos), value = TRUE)
  dados_comparativos[nir_cols] <- lapply(dados_comparativos[nir_cols], as.numeric)
  
  dados_grafico_comparativo <- dados_comparativos %>%
    mutate(Espectro_ID = paste(Specie, row_number(), sep = "_")) %>%
    tidyr::pivot_longer(cols = all_of(nir_cols), 
                        names_to = "Wavelength", 
                        values_to = "Absorbance") %>%
    mutate(Wavelength = as.numeric(sub("^X", "", Wavelength)))
  
  # Criar gráfico comparativo
  p_comparativo <- ggplot(dados_grafico_comparativo, aes(x = Wavelength, y = Absorbance, 
                                                         color = Specie, linetype = Metodo)) +
    geom_line(aes(group = interaction(Espectro_ID, Metodo)), linewidth = 0.5, alpha = 0.8) +
    scale_color_manual(values = paleta_cores, name = "Species") +
    scale_linetype_discrete(name = "Método de Pré-processamento") +
    labs(
      title = "Comparação de Métodos de Pré-processamento",
      x = expression("Wavenumber (cm"^{-1}*")"),
      y = "Absorbance"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    facet_wrap(~Metodo, ncol = 2, scales = "free_y")
  
  # Salvar gráfico comparativo
  ggsave("FIGURAS/Comparacao_Metodos_Preprocessamento.png",
         plot = p_comparativo, width = 12, height = 8, dpi = 300, units = "in")
}

cat("Processamento concluído! Verifique os gráficos na pasta FIGURAS.\n")




