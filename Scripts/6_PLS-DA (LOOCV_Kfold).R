# =============================================================================
# PLS-da (LOOCV/K-fold)
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
install.packages(c("dplyr", "caret", "MASS", "pls", "magrittr", "boot", "reshape2"))

library(dplyr)        # Para manipulação de dados
library(caret)        # Para machine learning
library(MASS)         # Para análises estatísticas
library(pls)          # Para regressão PLS
library(magrittr)     # Para usar o operador pipe %>% 
library(boot)         # Para métodos bootstrap
library(reshape2)     # Para transformar dados
# =============================================================================
# CRIE A TABELAS COM OS ARQUIVOS DE LEITURA
# =============================================================================

rm(list=ls()) #Limpar a lista de arquivos

setwd("C:/Users/nikso/OneDrive/Ms_métodos_NIR/ARTICLE_AND_ANALYSES/DIRÉTORIO_ANÁLISES")

#Carrega dado NIR
NIRdata = read.table ('TABELAS/DADOS MEDIOS_Aba.csv', header = TRUE, sep = ",", dec = ".", fileEncoding = "UTF-8", colClasses = "character")
NIRdata = read.table ('TABELAS/DADOS MEDIOS_Ada.csv', header = TRUE, sep = ",", dec = ".", fileEncoding = "UTF-8", colClasses = "character")
NIRdata = read.table ('TABELAS/DADOS MEDIOS_AbaAda.csv', header = TRUE, sep = ",", dec = ".", fileEncoding = "UTF-8", colClasses = "character")

NIRdata$Specie <- make.names(NIRdata$Specie)

dados = as.data.frame(NIRdata)

# Identificar as colunas espectrais
nir_cols <- 7:ncol(dados)

# Converter para numérico
dados[nir_cols] <- lapply(dados[nir_cols], function(x) as.numeric(as.character(x)))

# =============================================================================
# SÓ USAR SE EU FOR DIVIDIR OS DADOS EM TREINO/TESTE
# =============================================================================

#define um inicio para sorteios aleatorios em todos os processos (garante reprodutibilidade)
set.seed(123)

#carrega o pacote necessario
indice = createDataPartition(dados$Specie,
                             p = .7, #define an 70%/30% train/test split of the dataset
                             list = FALSE, 
                             times = 1)
treino = dados[indice,]
table(treino$Species)

teste = dados[-indice,]
table(teste$Species)

# =============================================================================
# PLS-DA USANDO A VALIDAÇÃO "LOOCV"
# =============================================================================

set.seed(123)
# Definindo controle para validação cruzada leave-one-out
pls_ctrl <- trainControl(
  method = "LOOCV",
  verboseIter = TRUE,
  classProbs = TRUE,
  returnData = TRUE,
  savePredictions = "final"
)

# Definindo um intervalo amplo para ncomp para permitir seleção automática
tuneLength = 20  # Testar até 20 componentes, o que deixa a escolha mais automática

pls_grid <- expand.grid(
  ncomp = seq(1, 20)  # Grid de valores para ncomp de 1 a 20
)

pls_loocv <- train(
  x = dados[, 7:ncol(dados)],
  y = dados[, 2],  
  method = "pls",
  trControl = pls_ctrl,
  tuneGrid = pls_grid,  # Utilizar a grade de parâmetros definida para ncomp
  tuneLength = tuneLength,  # Tentar ajustar com ncomp = 5 inicialmente
  verbose = TRUE,
  maxit = 4000  # Número máximo de iterações
)

# Criando a matriz de confusão para os dados de teste
conf_matrix_teste <- confusionMatrix(pls_loocv$pred$pred, pls_loocv$pred$obs)

# Exibindo os resultados da matriz de confusão
model_accuracy_teste <- conf_matrix_teste
print(model_accuracy_teste)

# Criando a matriz de confusão original
conf_matrix_teste <- confusionMatrix(pls_loocv$pred$pred, pls_loocv$pred$obs)

# Extraindo a TABELAS de confusão
conf_table <- conf_matrix_teste$table

# Transpondo a TABELAS de confusão
conf_table_transposed <- t(conf_table)

# Criando uma nova matriz de confusão com a TABELAS transposta
conf_matrix_teste_transposed <- confusionMatrix(conf_table_transposed)

# Exibindo os resultados da nova matriz de confusão
print(conf_matrix_teste_transposed)

# Criando um vetor com os 8 valores desejados para a nova coluna "BalAcc%"
BalAcc_values <- c(93, 90, 92, 92, 100, 100, 86, 89, 100, 93, 90, 99, 92)

# Adicionando os valores como uma nova coluna na matriz de confusão transposta
conf_table_transposed <- cbind(conf_table_transposed, "BalAcc%" = BalAcc_values)

# Criando um vetor com os nomes das classes
classes_values <- c("ADILPUL", "ANEHIR", "BLEOCC", "CYAMIC", "CYCMEN", "DICFLE", "DIPCRI", "LYGVEN", "MENSER", "MICGEM", "NIPCRA", "PLEAST", "SERTRI")

# Adicionando a coluna "Classes" como a primeira coluna
conf_table_transposed <- cbind(Classes = classes_values, conf_table_transposed)

# Convertendo a matriz em data frame e adicionando colunas
conf_table_transposed <- conf_table_transposed %>%
  as.data.frame(stringsAsFactors = FALSE, check.names = FALSE) %>%  # evita alteração de nomes
  cbind("BalAcc%" = BalAcc_values) %>%                              # adiciona BalAcc% fixo
  cbind(Classes = classes_values, .)                                # adiciona Classes como primeira coluna

# Garantindo a ordem das colunas
conf_table_transposed <- conf_table_transposed[, c("Classes", setdiff(colnames(conf_table_transposed), "Classes"))]

# Exibindo a matriz de confusão transposta atualizada sem aspas
print(conf_table_transposed, quote = FALSE)

# Exportando os resultados da matriz de confusão por classe
write.table(conf_table_transposed, file='TABELAS/Aba_matriz_loocv.csv', sep="\t", row.names=TRUE)
write.table(conf_table_transposed, file='TABELAS/Ada_matriz_loocv.csv', sep="\t", row.names=TRUE)
write.table(conf_table_transposed, file='TABELAS/AbaAda_matriz_loocv.csv', sep="\t", row.names=TRUE)

# =============================================================================
# PLS-DA USANDO A VALIDAÇÃO "K-FOLD"
# =============================================================================

rm(list=ls()) #Limpar a lista de arquivos

set.seed(123)
#Definindo controle para validação cruzada k-fold
pls_ctrl <- trainControl(
  method = "cv",  # Usando k-fold cross-validation
  number = 10,    # Número de folds (k)
  verboseIter = TRUE,
  classProbs = TRUE,
  returnData = TRUE,
  savePredictions = "final"
)

#Definindo um intervalo amplo para ncomp para permitir seleção automática
tuneLength = 20  # Testar até 20 componentes, o que deixa a escolha mais automática

#Treinando o modelo PLS-DA
pls_kfold <- train(
  x = dados[, 7:ncol(dados)],  # Selecione suas variáveis preditoras corretamente
  y = dados[, 2],               # Selecione sua variável de resposta corretamente
  method = "pls",
  trControl = pls_ctrl,
  tuneLength = tuneLength,  # Deixar o caret escolher o melhor ncomp automaticamente
  verbose = TRUE,
  maxit = 4000    # Número máximo de iterações
)

# Criando a matriz de confusão para os dados de teste
conf_matrix_teste <- confusionMatrix(pls_kfold$pred$pred, pls_kfold$pred$obs)

# Exibindo os resultados da matriz de confusão
model_accuracy_teste <- conf_matrix_teste
print(model_accuracy_teste)

# Criando a matriz de confusão original
conf_matrix_teste <- confusionMatrix(pls_kfold$pred$pred, pls_kfold$pred$obs)

# Extraindo a TABELAS de confusão
conf_table <- conf_matrix_teste$table

# Transpondo a TABELAS de confusão
conf_table_transposed <- t(conf_table)

# Criando uma nova matriz de confusão com a TABELAS transposta
conf_matrix_teste_transposed <- confusionMatrix(conf_table_transposed)

# Exibindo os resultados da nova matriz de confusão
print(conf_matrix_teste_transposed)

# Criando um vetor com os 8 valores desejados para a nova coluna "BalAcc%"
BalAcc_values <- c(93, 95, 93, 95, 100, 100, 89, 95, 100, 93, 90, 99, 92)

# Adicionando os valores como uma nova coluna na matriz de confusão transposta
conf_table_transposed <- cbind(conf_table_transposed, "BalAcc%" = BalAcc_values)

# Criando um vetor com os nomes das classes
classes_values <- c("ADILPUL", "ANEHIR", "BLEOCC", "CYAMIC", "CYCMEN", "DICFLE", "DIPCRI", "LYGVEN", "MENSER", "MICGEM", "NIPCRA", "PLEAST", "SERTRI")

# Adicionando a coluna "Classes" como a primeira coluna
conf_table_transposed <- cbind(Classes = classes_values, conf_table_transposed)

# Convertendo a matriz em data frame e adicionando colunas
conf_table_transposed <- conf_table_transposed %>%
  as.data.frame(stringsAsFactors = FALSE, check.names = FALSE) %>%  # evita alteração de nomes
  cbind("BalAcc%" = BalAcc_values) %>%                              # adiciona BalAcc% fixo
  cbind(Classes = classes_values, .)                                # adiciona Classes como primeira coluna

# Garantindo a ordem das colunas
conf_table_transposed <- conf_table_transposed[, c("Classes", setdiff(colnames(conf_table_transposed), "Classes"))]

# Exibindo a matriz de confusão transposta atualizada sem aspas
print(conf_table_transposed, quote = FALSE)

# Exportando os resultados da matriz de confusão por classe
write.table(conf_table_transposed, file='TABELAS/Aba_matriz_K-fold.csv', sep="\t", row.names=TRUE)
write.table(conf_table_transposed, file='TABELAS/Ada_matriz_K-fold.csv', sep="\t", row.names=TRUE)
write.table(conf_table_transposed, file='TABELAS/AbaAda_matriz_K-fold.csv', sep="\t", row.names=TRUE)

#----------------------------------------------------------------------------------------------------------------------------
rm(list=ls()) #Limpar a lista de arquivos

# Para plotar graficos

# LOOCV
tab = read.table('TABELAS/Aba_matriz_loocv.csv',sep = "\t", stringsAsFactors = TRUE, header = TRUE, check.names = FALSE)
tab = read.table('TABELAS/Ada_matriz_loocv.csv',sep = "\t", stringsAsFactors = TRUE, header = TRUE, check.names = FALSE)
tab = read.table('TABELAS/AbaAda_matriz_loocv.csv',sep = "\t", stringsAsFactors = TRUE, header = TRUE, check.names = FALSE)

# K-FOLD
tab = read.table('TABELAS/Aba_matriz_K-fold.csv',sep = "\t", stringsAsFactors = TRUE, header = TRUE, check.names = FALSE)
tab <- read.table('TABELAS/Ada_matriz_K-fold.csv', sep = "\t", stringsAsFactors = TRUE, header = TRUE, check.names = FALSE)
tab = read.table('TABELAS/AbaAda_matriz_K-fold.csv',sep = "\t", stringsAsFactors = TRUE, header = TRUE, check.names = FALSE)

#Transformar a TABELAS em um formato longo
data_long <- melt(tab, id.vars = "Classes")

#Renomear as colunas
colnames(data_long) <- c("Classe_predita", "Classe_real", "Valor")

#Criar a matriz de confusão
confusion_matrix <- table(data_long$Classe_predita, data_long$Classe_real)

# ==============================================================================
# Opção 1: Matriz de confusão convencional
# ==============================================================================

# Plot com paleta azul escuro mas fundo claro
ggplot() +
  geom_tile(data = filter(data_long, Classe_real != "BalAcc%" & Valor != 0), 
            aes(x = Classe_real, y = Classe_predita, fill = Valor), 
            color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("#2C5282", "#3182CE", "#90CDF4"),  # Gradiente azul escuro
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = filter(data_long, Classe_real == "BalAcc%"), 
            aes(x = Classe_real, y = Classe_predita, fill = "highlight"), 
            color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("highlight" = "#1A365D"), guide = "none") +
  geom_text(data = filter(data_long, Valor != 0), 
            aes(x = Classe_real, y = Classe_predita, 
                label = ifelse(Valor != 0, Valor, ""), 
                color = "white",  # Todos os textos em branco
                fontface = ifelse(Classe_real == "BalAcc%", "bold", "plain")), 
            size = 5) +
  scale_color_manual(values = c("white"), guide = "none") +  # Apenas cor branca
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.grid = element_blank()
  ) +
  labs(x = "", y = "")

# Salvar o gráfico como um arquivo PNG

# LOOCV
ggsave("Figuras/Aba_matrizconfusao_fertil_loocv.png", width = 7, height = 7, dpi = 300)
ggsave("Figuras/Ada_matrizconfusao_esteril_loocv.png", width = 7, height = 7, dpi = 300)
ggsave("Figuras/AbaAda_matrizconfusao_combinada_loocv.png", width = 7, height = 7, dpi = 300)

# K-FOLD
ggsave("Figuras/Aba_matrizconfusao_fertil_K-fold.png", width = 7, height = 7, dpi = 300)
ggsave("Figuras/Ada_matrizconfusao_esteril_K-fold.png", width = 7, height = 7, dpi = 300)
ggsave("Figuras/AbaAda_matrizconfusao_combinada_K-fold.png", width = 7, height = 7, dpi = 300)
