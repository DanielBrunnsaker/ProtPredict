library(clusterProfiler)
library(AnnotationDbi)
library(org.Sc.sgd.db)
library(patchwork)

#load the gene lists?
background = read.csv('/Users/danbru/Desktop/MetPredict/ProteinList.tsv', sep = '\t') #all 2292 proteins
background = read.csv('/Users/danbru/Desktop/MetPredict/predictableproteins.tsv', sep = '\t') #all predictable proterins


proteins = read.csv('/Users/danbru/Desktop/MetPredict/AAPredictability.tsv', sep = '\t')
glutamine = proteins$glutamine[proteins$glutamine != '']
alanine = proteins$alanine[proteins$alanine != '']
glycine = proteins$glycine[proteins$glycine != '']
proline = proteins$proline[proteins$proline != '']
arginine = proteins$arginine[proteins$arginine != '']
serine = proteins$serine[proteins$serine != '']
tryptophan = proteins$tryptophan[proteins$tryptophan != '']


glutamine_entrez = select(org.Sc.sgd.db, glutamine, "ENTREZID", "UNIPROT")
alanine_entrez = select(org.Sc.sgd.db, alanine, "ENTREZID", "UNIPROT")
glycine_entrez = select(org.Sc.sgd.db, glycine, "ENTREZID", "UNIPROT")
proline_entrez = select(org.Sc.sgd.db, proline, "ENTREZID", "UNIPROT")
arginine_entrez = select(org.Sc.sgd.db, arginine, "ENTREZID", "UNIPROT")
serine_entrez = select(org.Sc.sgd.db, serine, "ENTREZID", "UNIPROT")
tryptophan_entrez = select(org.Sc.sgd.db, tryptophan, "ENTREZID", "UNIPROT")


background_entrez = select(org.Sc.sgd.db, background$X0, "ENTREZID", "UNIPROT")

#glutamine
results_BP = enrichGO(glutamine_entrez$ENTREZID,'org.Sc.sgd.db', ont="BP", universe = background_entrez$ENTREZID)
results_MF = enrichGO(glutamine_entrez$ENTREZID,'org.Sc.sgd.db', ont="MF", universe = background_entrez$ENTREZID)
results_CC = enrichGO(glutamine_entrez$ENTREZID,'org.Sc.sgd.db', ont="CC", universe = background_entrez$ENTREZID)

bp = dotplot(results_BP)
mf = dotplot(results_MF)
cc = dotplot(results_CC)

bp+cc

#alanine, no enrichment, will bug out
results_BP = enrichGO(alanine_entrez$ENTREZID,'org.Sc.sgd.db', ont="BP", universe = background_entrez$ENTREZID)
results_MF = enrichGO(alanine_entrez$ENTREZID,'org.Sc.sgd.db', ont="MF", universe = background_entrez$ENTREZID)
results_CC = enrichGO(alanine_entrez$ENTREZID,'org.Sc.sgd.db', ont="CC", universe = background_entrez$ENTREZID)

p1 = dotplot(results_BP)
p2 = dotplot(results_MF)
p3 = dotplot(results_CC)

p1+p2+p3

#glycine
results_BP = enrichGO(glycine_entrez$ENTREZID,'org.Sc.sgd.db', ont="BP", universe = background_entrez$ENTREZID)
results_MF = enrichGO(glycine_entrez$ENTREZID,'org.Sc.sgd.db', ont="MF", universe = background_entrez$ENTREZID)
results_CC = enrichGO(glycine_entrez$ENTREZID,'org.Sc.sgd.db', ont="CC", universe = background_entrez$ENTREZID)

p1 = dotplot(results_BP)
p2 = dotplot(results_MF)
p3 = dotplot(results_CC)

p1+p2+p3

#proline
results_BP = enrichGO(proline_entrez$ENTREZID,'org.Sc.sgd.db', ont="BP", universe = background_entrez$ENTREZID)
results_MF = enrichGO(proline_entrez$ENTREZID,'org.Sc.sgd.db', ont="MF", universe = background_entrez$ENTREZID)
results_CC = enrichGO(proline_entrez$ENTREZID,'org.Sc.sgd.db', ont="CC", universe = background_entrez$ENTREZID)

p1 = dotplot(results_BP)
p2 = dotplot(results_MF)
p3 = dotplot(results_CC)

p1+p2+p3

#arginine, null.
results_BP = enrichGO(arginine_entrez$ENTREZID,'org.Sc.sgd.db', ont="BP", universe = background_entrez$ENTREZID)
results_MF = enrichGO(arginine_entrez$ENTREZID,'org.Sc.sgd.db', ont="MF", universe = background_entrez$ENTREZID)
results_CC = enrichGO(arginine_entrez$ENTREZID,'org.Sc.sgd.db', ont="CC", universe = background_entrez$ENTREZID)

p1 = dotplot(results_BP)
p2 = dotplot(results_MF)
p3 = dotplot(results_CC)

p1+p2+p3

#serine
results_BP = enrichGO(serine_entrez$ENTREZID,'org.Sc.sgd.db', ont="BP", universe = background_entrez$ENTREZID)
results_MF = enrichGO(serine_entrez$ENTREZID,'org.Sc.sgd.db', ont="MF", universe = background_entrez$ENTREZID)
results_CC = enrichGO(serine_entrez$ENTREZID,'org.Sc.sgd.db', ont="CC", universe = background_entrez$ENTREZID)

p1 = dotplot(results_BP)
p2 = dotplot(results_MF)
p3 = dotplot(results_CC)

p1+p2+p3


#tryptophan, null
results_BP = enrichGO(tryptophan_entrez$ENTREZID,'org.Sc.sgd.db', ont="BP", universe = background_entrez$ENTREZID)
results_MF = enrichGO(tryptophan_entrez$ENTREZID,'org.Sc.sgd.db', ont="MF", universe = background_entrez$ENTREZID)
results_CC = enrichGO(tryptophan_entrez$ENTREZID,'org.Sc.sgd.db', ont="CC", universe = background_entrez$ENTREZID)

p1 = dotplot(results_BP)
p2 = dotplot(results_MF)
p3 = dotplot(results_CC)

p1+p2+p3


#let's roll with glutamine, glycine  and proline




