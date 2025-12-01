# =============================================================================
# CRIAR TABELAS DE MODELOS COM MÉDIAS (Aba, Ada e AbaAda)
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
install.packages("dplyr")

library(dplyr) # Manipulação de dados (filtros, seleção, agrupamento)
# =============================================================================
# FILTRANDO DADOS PARA CRIAR TABELAS
# =============================================================================

# Limpa o ambiente e fecha todas as visualizações
rm(list = ls()) 

# Diretorio das tabelas
setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES/")

# Importar a tabela dados NIR
Matriz_NIR <- read.table(file='TABELAS/Matriz_NIR_snv.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

# Exluir numeros dedepois de ABA ou ADA
Matriz_NIR <- Matriz_NIR %>%
  filter(substr(Abaxial_or_Adaxial, 1, 3) %in% c("ABA", "ADA")) %>%
  mutate(Abaxial_or_Adaxial = substr(Abaxial_or_Adaxial, 1, 3))  # Remove números

# Definir coluna e classes para filtrar
quecoluna <- "Abaxial_or_Adaxial"  # Nome da coluna
queclasse1 <- "ABA"                # Filtro 1 (3 primeiras letras)
queclasse2 <- "ADA"                # Filtro 2 (3 primeiras letras)

# cria um vetor logico, ou seja TRUE/FALSE para cada linha que corresponde ou não ao filtro
vl <- Matriz_NIR[,quecoluna]==queclasse1
v2 <- Matriz_NIR[,quecoluna]==queclasse2

# Filtrar os dados
dados.abaxialbrutos <- Matriz_NIR[vl, ]
dados.adaxialbrutos <- Matriz_NIR[v2, ]

# Combinar os dados filtrados em uma única tabela
dados.combinados <- rbind(dados.abaxialbrutos, dados.adaxialbrutos)

# salva o dado filtrado
write.table(dados.abaxialbrutos, file='TABELAS/DADOS BRUTOS_abaxial.csv',row.names = FALSE, sep = ";", quote = FALSE, dec = ".", fileEncoding = "UTF-8")              
write.table(dados.adaxialbrutos, file='TABELAS/DADOS BRUTOS_adaxial.csv',row.names = FALSE, sep = ";", quote = FALSE, dec = ".", fileEncoding = "UTF-8") 
write.table(dados.combinados, file='TABELAS/DADOS BRUTOS_AbaAda.csv',row.names = FALSE, sep = ";", quote = FALSE, dec = ".", fileEncoding = "UTF-8") 

# =============================================================================
# FILTRANDO DADOS PARA CRIAR TABELAS COM DADOS MÉDIOS
# =============================================================================

rm(list=ls()) #Limpar a lista de arquivos

# Importar a tabela com dados brutos 3001
dados_Aba <- read.table(file='TABELAS/DADOS BRUTOS_abaxial.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")
dados_Ada <- read.table(file='TABELAS/DADOS BRUTOS_adaxial.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")
dados_AbaAda <- read.table(file='TABELAS/DADOS BRUTOS_AbaAda.csv',header = TRUE, sep = ";", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

colspecie <- "Specie"
colcollector <- "Collector_and_number"
colFE <- "Fertile_or_Sterile"
colDM <- "Dimorphic_or_Monomorphic"
colFnumb <- "Frond_Number"
colABAD <- "Abaxial_or_Adaxial"

# 3. Lista com os dataframes e seus nomes de saída
lista_dados <- list(
  "Aba" = dados_Aba,
  "Ada" = dados_Ada,
  "AbaAda" = dados_AbaAda)

# 4. Processamento automatizado para cada dataframe
for (nome in names(lista_dados)) {
  dados <- lista_dados[[nome]]
  
  # Selecionar colunas NIR (X...)
  colunas_nir <- grep("^X", colnames(dados), value = TRUE)
  
  # Converter para numérico e calcular média por coletor
  dados_nir <- as.data.frame(lapply(dados[, colunas_nir], as.numeric))
  medias_nir <- aggregate(dados_nir, by = list(Coletor = dados[[colcollector]]), FUN = mean, na.rm = TRUE)
  colnames(medias_nir)[1] <- colcollector
  
  # Extrair metadados únicos
  metadados <- dados[, c(colspecie, colcollector, colFE, colDM, colFnumb, colABAD)]
  metadados <- metadados[!duplicated(metadados[[colcollector]]), ]
  
  # Juntar tudo
  dados_processados <- merge(metadados, medias_nir, by = colcollector)
  
  # Salvar em arquivo
  nome_arquivo <- paste0("TABELAS/DADOS MEDIOS_", nome, ".csv")
  write.csv(dados_processados, nome_arquivo, row.names = FALSE)
  
  cat("Arquivo salvo:", nome_arquivo, "\n")
}

# Após o loop, carregue os arquivos salvos 3001:
View(read.csv("TABELAS/DADOS MEDIOS_Aba.csv"))
View(read.csv("TABELAS/DADOS MEDIOS_Ada.csv"))
View(read.csv("TABELAS/DADOS MEDIOS_AbaAda.csv"))
