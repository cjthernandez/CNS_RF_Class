library(VennDiagram)

# These RF models need to be generated
load(file.path("Results","rf.pred_Capper.RData"))
imp_cap <- rf.pred$importance
load(file.path("Results","rf.pred_32k.RData"))
imp_32k <- rf.pred$importance
load(file.path("Results","rf.pred_500trees.RData"))
imp_500 <- rf.pred$importance
load(file.path("Results","rf.pred_32k500trees.RData"))
imp_32k500 <- rf.pred$importance


probes <- c(rownames(imp_cap), rownames(imp_32k), rownames(imp_500), rownames(imp_32k500))
imps <- c(imp_cap[,92], imp_32k[,92], imp_500[,92], imp_32k500[,92]) # Permutation importance
mods <- c(rep("Capper", nrow(imp_cap)), rep("32k", nrow(imp_32k)), rep("500t", nrow(imp_500)), rep("32k500t", nrow(imp_32k500)))

df_imps <- data.frame(probes,imps,mods)
df_ii <- df_imps %>% group_by(probes) %>% mutate(shared=n()) %>% ungroup()

# Permutation importance of all probes in each model vs probes shared between all models
ggplot(df_ii, aes(x=imps)) + 
  geom_histogram(data=df_ii, aes(y=..count../sum(..count..)), fill="blue",alpha=0.2) +
  geom_histogram(data=subset(df_ii, shared==4), aes(y=..count../sum(..count..)), fill="red",alpha=0.2) +
  facet_wrap(~mods)


venn.diagram(list(Capper=rownames(imp_cap), Capper32k=rownames(imp_32k), Capper500t=rownames(imp_500), Capper32k500t=rownames(imp_32k500)),
             filename="Rplots/venn_diagram_Capper_probes.tiff")