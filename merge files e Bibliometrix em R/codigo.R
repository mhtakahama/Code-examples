#install.packages("bibliometrix")

library(bibliometrix) # carregar o pacote na memória para utilizá-lo agora.

# Obtém o diretório do script R atual
diretorio_atual <- getwd()
diretorio_atual2 <- setwd(diretorio_atual)

# Constrói o caminho completo para o arquivo savedrecs.bib
caminho_arquivo1 <- file.path(diretorio_atual, "savedrecs.bib")
# Constrói o caminho completo para o arquivo savedrecs.bib
caminho_arquivo2 <- file.path(diretorio_atual, "scopus.bib")

# SCOPUS: Converter os dados   para o padrão do bibliometrix
A <- convert2df(caminho_arquivo1, dbsource = "isi", format = "bibtex")
B <- convert2df(caminho_arquivo2, dbsource = "scopus", format = "bibtex")
# C <- convert2df("c:/bib/citation-export.bib", dbsource = "isi", format = "bibtex")

M <- mergeDbSources(A, B, remove.duplicated = TRUE)
# M <- mergeDbSources(A, B, C, remove.duplicated = TRUE) 

#### Cria um arquivo.csv para importar para o Excel  
P<- M[,c("AU","TI","AB","DE","ID","SO","TC","PY","LA","DT","DI")]  # Cria lista na ordem desejada

caminho_arquivo3 <- file.path(diretorio_atual, "artigos_conjuntos.csv")
write.table(M, caminho_arquivo3, sep=";", row.names=FALSE) # Para gerar com separador ";", sem necessidade de ajustar o CSV, quanto à primeira coluna

# Criar Dados conjuntos 
write.table(M, "c:/bib/artigos_conjuntos.csv", sep=";", row.names=FALSE) # Para usar no dados completos de fontes no biblioshiny()


resultados <- biblioAnalysis(M) 
Resumo <- summary(object = resultados, k = 10)
biblioshiny()