library(topGO)

# Construct the file paths using the provided parameter
interesting_genes_path <- commandArgs(trailingOnly = TRUE)[1] #geneList.txt
all_genes_path <- commandArgs(trailingOnly = TRUE)[2] #all_genes.txt
topGO_formatted_path <- commandArgs(trailingOnly = TRUE)[3] #GO_formatted.txt
output_path <- commandArgs(trailingOnly = TRUE)[4] #output.txt

# Read gene lists
interesting_genes <- readLines(interesting_genes_path)
genes_all <- readLines(all_genes_path)

all_genes <- unique(c(genes_all, interesting_genes))

# Create a factor with levels indicating the source file
gene_factor <- factor(all_genes, levels = all_genes)

# Assign levels to indicate the source file (0 for not interesting, 1 for interesting)
gene_levels <- ifelse(all_genes %in% interesting_genes, 1, 0)

# Assign levels to the factor
levels(gene_factor) <- gene_levels

# Assign gene names to the factor
names(gene_factor) <- all_genes

geneID2GO <- readMappings(file = topGO_formatted_path)

BP_GOdata <- new("topGOdata", ontology = "MF", allGenes = gene_factor, annot = annFUN.gene2GO, gene2GO = geneID2GO)

test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(BP_GOdata, test.stat)
resultData <- geneData(resultFisher)

allRes <- GenTable(BP_GOdata, classic = resultFisher, #KS = resultKS, weight = resultWeight,
                  orderBy = "weight", ranksOf = "classic", topNodes = 21)

data <- data.frame(GO.ID = allRes$GO.ID,
                   Term = allRes$Term,
                   Annotated = allRes$Annotated,
                   Significant = allRes$Significant,
                   Expected = allRes$Expected,
                   Fisher = allRes$classic)

# Save the final table
write.table(data, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE, )