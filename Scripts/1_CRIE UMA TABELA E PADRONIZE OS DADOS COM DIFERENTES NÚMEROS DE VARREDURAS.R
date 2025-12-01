# =============================================================================
# CRIE UMA TABELA E PADRONIZE OS DADOS COM DIFERENTES NÚMEROS DE VARREDURAS
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
install.packages("readr", "dplyr")

library(readr)  # Leitura eficiente de dados (CSV, TXT)
library(dplyr)  # Manipulação de dados (filtros, seleção, agrupamento)
# =============================================================================
# CRIE A TABELA COM OS ARQUIVOS DE LEITURA
# =============================================================================

#Diretório de onde estão os arquivos .csv dos espectros
folder_path <- "C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/ESPECTROS_19072025"

files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

result_list <- list()

for (file_path in files) {
  df <- read.csv(file_path, sep = ";", stringsAsFactors = FALSE, skip = 1)
  
  if(ncol(df) >= 2) {
    # Obter os comprimentos de onda e os valores de absorbância
    cm1 <- df[, 1]
    absorb <- df[, 2]
    
    # Corrige possíveis vírgulas para pontos
    absorb <- gsub(",", ".", absorb)
    absorb <- as.numeric(absorb)
    
    # Nomes das colunas com base nos valores de cm-1
    nomes_colunas <- as.character(cm1)
    
    # Processa nome do arquivo
    file_name <- tools::file_path_sans_ext(basename(file_path))
    file_name <- gsub("A$", "", file_name)
    file_name <- gsub("\\.$", "", file_name)
    file_parts <- unlist(strsplit(file_name, "_"))
    
    if (length(file_parts) >= 4) {
      fertile_letter <- substr(file_parts[3], 1, 1)
      dimorphic_letter <- substr(file_parts[3], 2, 2)
      frond_number <- substr(file_parts[3], 3, nchar(file_parts[3]))
      
      # Cria metadados
      temp_df <- data.frame(
        Specie = file_parts[1],
        Collector_and_number = file_parts[2],
        Fertile_or_Sterile = fertile_letter,
        Dimorphic_or_Monomorphic = dimorphic_letter,
        Frond_Number = frond_number,
        Abaxial_or_Adaxial = file_parts[4],
        stringsAsFactors = FALSE
      )
      
      # Cria a linha de absorbância com nomes dos comprimentos de onda
      temp_abs <- as.data.frame(t(absorb))
      colnames(temp_abs) <- nomes_colunas
      
      # Junta tudo
      temp_df <- cbind(temp_df, temp_abs)
      
      # Nome do grupo conforme quantidade de pontos espectrais
      nome_grupo <- paste0("absorbancia_", length(absorb), "_valores")
      
      if (!nome_grupo %in% names(result_list)) {
        result_list[[nome_grupo]] <- temp_df
      } else {
        result_list[[nome_grupo]] <- rbind(result_list[[nome_grupo]], temp_df)
      }
    }
  }
}

# Salvar arquivos separados
for (nome in names(result_list)) {
  df <- result_list[[nome]]
  
  output_path <- file.path("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/TABELAS", paste0("TAB_", nome, ".csv"))
  
  write.table(df,
              file = output_path,
              row.names = FALSE,
              sep = ";",
              quote = FALSE,
              dec = ".",
              fileEncoding = "UTF-8")
}

cat("✅ Tabelas geradas com colunas nomeadas por valores cm-1.\n")

# =============================================================================
# AJUSTE AS TABELAS DEVIDO AOS DIFERENTES VALORES DE VARREDURA 
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/TABELAS")

# Leitura do arquivo CSV
df_lido3001 <- read.table("TAB_absorbancia_3001_valores.csv", header = TRUE, 
                          sep = ";", dec = ".", fileEncoding = "UTF-8", 
                          colClasses = "character")

df_lido6001 <- read.table("TAB_absorbancia_6001_valores.csv", header = TRUE, 
                          sep = ";", dec = ".", fileEncoding = "UTF-8", 
                          colClasses = "character")

# Corrigir nomes das colunas da 7 em diante: manter o "X", e remover ".0", ".00"
colnames(df_lido3001) <- gsub("\\.0+$", "", colnames(df_lido3001))
colnames(df_lido6001) <- gsub("\\.0+$", "", colnames(df_lido6001))

# Definir as colunas de referência
colunas_referencia <- colnames(df_lido3001)
colunas_fixas <- colnames(df_lido3001)[1:6]

# Função para ajustar os dataframes
ajustar_dataframe <- function(df, ref_cols, fixed_cols) {
  # Filtrar colunas variáveis (a partir da 7ª)
  colunas_variaveis <- colnames(df)[7:ncol(df)]
  colunas_variaveis_para_manter <- colunas_variaveis[colunas_variaveis %in% ref_cols]
  
  # Selecionar colunas para manter
  colunas_para_manter <- c(fixed_cols, colunas_variaveis_para_manter)
  df_ajustado <- df[, colunas_para_manter, drop = FALSE]
  
  # Ordenar colunas variáveis conforme referência
  ordem_variaveis <- ref_cols[ref_cols %in% colunas_variaveis_para_manter]
  df_ajustado <- df_ajustado[, c(fixed_cols, ordem_variaveis)]
  
  return(df_ajustado)
}

# Aplicar a função
df_lido6001 <- ajustar_dataframe(df_lido6001, colunas_referencia, colunas_fixas)

# Verificar se os nomes das colunas são iguais
identical(colnames(df_lido3001), colnames(df_lido6001))

# Unir todos os dataframes em um único dataframe
df_final <- rbind(df_lido3001, df_lido6001)  

#Antes de salvar eu vou recortar quero so essa região (4000 a 7000 cm-1)
df_final <- df_final[, -c(7:1506)]

# Salvar arquivo
write.table(df_final,                             # Objeto dataframe a ser salvo
            file = "Matriz_NIR.csv",              # Nome do arquivo de saída
            row.names = FALSE,                    # Não incluir números de linha
            sep = ";",                            # Separador de colunas
            quote = FALSE,                        # Não usar aspas nos valores
            dec = ".",                            # Ponto como separador decimal
            fileEncoding = "UTF-8")               # Codificação do arquivo

# =============================================================================
# CALCULAR QUANTOS ESPÉCIMES TEM POR ESPÉCIE
# =============================================================================

# Contar leituras únicas por espécie (quantos espécimes por espécie tem?)
contagem <- Matriz_NIR %>%
  distinct(Specie, Collector_and_number) %>%  # Remove duplicatas de combinações
  count(Specie, name = "Leituras_Unicas")     # Conta ocorrências por espécie

# Mostrar o resultado
print(contagem)
