# msig_to_setrankDb = function(msigdb, orga, dbname){
#   # msigdb : database R object
#   # orga : organism name (HUMAN|MOUSE)
#   # dbname : Msig database name (C2|H)
#   i = 1
#   df = data.frame(geneID = factor(), termID = factor(), termName = factor(), dbName = factor(), description = factor())
#   for (name in(names(msigdb))){
#     for (gene in msigdb[[name]]){
#       df = rbind(df, data.frame(geneID = factor(gene), termID = factor(paste(orga,"MSigDB",dbname,i,sep = "_")), termName = factor(name), dbName = factor("MSigDB"), description = factor("")))
#     }
#     i = i+1
#   }
#   return(df)
# }
#
# # #MOUSE
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_H_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c2_v5p2.rdata"))
# #
# #
# # #HUMAN
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata"))
# load(url("http://bioinf.wehi.edu.au/software/MSigDB/human_c2_v5p2.rdata"))
#
#
# Mm.c2.kegg.subset = c(Mm.c2[grep("KEGG",attributes(Mm.c2)$names)])
# Mm.c2.reactome.subset = c(Mm.c2[grep("REACTOME", attributes(Mm.c2)$names)])
#
# Hs.c2.kegg.subset = c(Hs.c2[grep("KEGG",attributes(Hs.c2)$names)])
# Hs.c2.reactome.subset = c(Hs.c2[grep("REACTOME", attributes(Hs.c2)$names)])
#
# mouse_H = msig_to_setrankDb(Mm.H, "MOUSE", "H")
# mouse_C2_kegg = msig_to_setrankDb(Mm.c2.kegg.subset,"MOUSE","C2Kegg")
# mouse_C2_reactome = msig_to_setrankDb(Mm.c2.reactome.subset,"MOUSE","C2Reactome")
#
#
#
# human_H = msig_to_setrankDb(Hs.H,"HUMAN","H")
# human_C2_kegg = msig_to_setrankDb(Hs.c2.kegg.subset,"HUMAN","C2Kegg")
# human_C2_reactome = msig_to_setrankDb(Hs.c2.reactome.subset,"HUMAN","C2Reactome")
#
# MusCollectionH = buildSetCollection(mouse_H, referenceSet = allMusGenes, maxSetSize = 1450)
# MusCollectionC2Kegg = buildSetCollection(mouse_C2_kegg, referenceSet = allMusGenes, maxSetSize = 1450)
# MusCollectionC2Reactome = buildSetCollection(mouse_C2_reactome, referenceSet = allMusGenes, maxSetSize = 1450)
#
# HomoCollectionH = buildSetCollection(human_H, referenceSet = allHomoGenes, maxSetSize = 1450)
# HomoCollectionC2Kegg = buildSetCollection(human_C2_kegg, referenceSet = allHomoGenes, maxSetSize = 1450)
# HomoCollectionC2Reactome = buildSetCollection(human_C2_reactome, referenceSet = allHomoGenes, maxSetSize = 1450)
#
#
# save(MusCollectionH, file = "setRankMusCollection_H.RData")
# save(MusCollectionC2Kegg, file = "setRankMusCollection_C2Kegg.RData")
# save(MusCollectionC2Reactome, file = "setRankMusCollection_C2Reactome.RData")
#
# save(HomoCollectionH, file = "setRankHomoCollection_H.RData")
# save(HomoCollectionC2Kegg, file = "setRankHomoCollection_C2Kegg.RData")
# save(HomoCollectionC2Reactome, file = "setRankHomoCollection_C2Reactome.RData")
#
# # load("~/GitLab/DEA-DifferentialExpressionAnalysis/setRankHomoCollection_H.RData")
# # load("~/GitLab/DEA-DifferentialExpressionAnalysis/setRankHomoCollection_C2Kegg.RData")
# # load("~/GitLab/DEA-DifferentialExpressionAnalysis/setRankHomoCollection_C2Reactome.RData")
# #
# # load("~/GitLab/DEA-DifferentialExpressionAnalysis/setRankMusCollection_H.RData")
# # load("~/GitLab/DEA-DifferentialExpressionAnalysis/setRankMusCollection_C2Kegg.RData")
# # load("~/GitLab/DEA-DifferentialExpressionAnalysis/setRankMusCollection_C2Reactome.RData")
#
# # too slow when installin package - just put the .rdata files in the /data directory
# # devtools::use_data(MusCollectionH, MusCollectionC2Kegg, MusCollectionC2Reactome, HomoCollectionH, HomoCollectionC2Kegg, HomoCollectionC2Reactome, internal = TRUE, compress = "gzip", overwrite = TRUE)
