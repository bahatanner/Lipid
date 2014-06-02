# get lipid species from lipid maps database (.sdf file)

setwd("D:/Dropbox/Eigene Dateien/Desktop/LMSDFDownload18Mar14")

lipidmaps_sdf<-read.delim("LMSDFDownload18Mar14FinalAll.sdf")

# get lipid species with KEGG_ID

KEGG_ID<-subset(lipidmaps_sdf, grepl("KEGG_ID", lipidmaps_sdf$LMFA00000001))
KEGG_ID=as.numeric(row.names(KEGG_ID))+1
KEGG_ID<-lipidmaps_sdf[KEGG_ID,]

# get lipid species with CHEBI_ID

CHEBI_ID<-subset(lipidmaps_sdf, grepl("CHEBI_ID", lipidmaps_sdf$LMFA00000001))
CHEBI_ID=as.numeric(row.names(CHEBI_ID))+1
CHEBI_ID<-paste("CHEBI:", lipidmaps_sdf[CHEBI_ID,], sep="")

# map lipid KEGG_ID and lipid CHEBI_ID to total yeast compounds in metabolic model

setwd("D:/Dropbox/Eigene Dateien/Desktop/Lipid/yeast_7_11")

# yeast compounds and enzymes

yeast_compounds_enzymes<-read.delim("species_par_lipid711.tsv", header=TRUE, sep="\t")
yeast_compounds<-subset(yeast_compounds_enzymes, grepl("s_", yeast_compounds_enzymes$SpeciesID))
yeast_KEGG<-subset(yeast_compounds, grepl("kegg.compound", yeast_compounds$Annotation))
yeast_CHEBI<-subset(yeast_compounds, grepl("CHEBI", yeast_compounds$Annotation))

# lipid KEGG_ID in yeast metabolic network

lipid_species_KEGG_which<-NULL
for(i in KEGG_ID){
  lipid_species_KEGG_which<-append(lipid_species_KEGG_which, which(yeast_KEGG$Annotation==paste("http://identifiers.org/kegg.compound/",i,sep="")))
}
lipid_species_KEGG<-yeast_KEGG[lipid_species_KEGG_which,]

write.table(lipid_species_KEGG, file="lipid_species_KEGG.txt", sep="\t")

# lipid CHEBI_ID in yeast metabolic network

lipid_species_CHEBI_which<-NULL
for(i in CHEBI_ID){
  lipid_species_CHEBI_which<-append(lipid_species_CHEBI_which, which(yeast_CHEBI$Annotation==paste("http://identifiers.org/chebi/",i,sep="")))
}
lipid_species_CHEBI<-yeast_CHEBI[lipid_species_CHEBI_which,]

write.table(lipid_species_CHEBI, file="lipid_species_CHEBI.txt", sep="\t")

# produce a lipid masterfile

lipid_species<-rbind(lipid_species_CHEBI, lipid_species_KEGG)
lipid_species<-lipid_species[!duplicated(lipid_species$SpeciesID), ]

write.table(lipid_species, file="lipid_species.txt", sep="\t")

# get lipid species by list

setwd("D:/Dropbox/Eigene Dateien/Desktop/Lipid/yeast_7_11")

rxn_yeast<-read.delim("rxn_lipid711.tsv", header=TRUE, sep="\t")

lipid_list<-unique(c("HETE", "mevalonic", "carnitine", "mevalonate", "epoxysqualene", "squalene", "lauroyl", "palmitoyl", "decanoyl", "hexacosanoyl", "CoA", "myo-inositol", "inositol", "glycerol", "acylglycero", "geranyl", "farnesyl", "sphinganine", "sphingosine", "sphingolipid", "phosphoserine", "acetyl", "CDP-choline", "CDP-ethanolamine", "ceramide", "choline", "serine", "ethanolamine", "prenyl", "leukotriene", "lipid", "fatty", "diglyceride", "lipoamide", "lipoyl", "dihydroxyacetone", "sterol", "ergosta", "farnesyl", "hexadecanal", "pentenyl", "laurate", "malonyl", "myristate", "myristoyl", "octanoate", "octanoyl", "oleate", "oleoyl", "palmitate", "propionyl", "stearate", "stearoyl", "tetracosanoyl", "enoyl", "triglyceride", "glyceride", "pantothenate", "oxopentanoate", "acyl", "palmitoleate", "ergosta", "malonyl", "coenzyme A", "tetradecanoyl", "icosanoyl", "docosanoyl", "hexacosanoyl", "palmitoleoyl", "oleoyl", "oleate", "myristate", "laurate", "lauroyl", "palmitoleoyl", "lignoceric", "icosanoyl", "docosanoyl", "dihydroxyacetone", "steryl", "ester", "phosphatidyl", "phosphatidate", "cerotic"))

lipid_species<-rxn_yeast[which(grepl(paste(lipid_list,collapse="|"), rxn_yeast$MetName)),]
lipid_species<-lipid_species[!duplicated(lipid_species$MetName), ]

write.table(lipid_species, file="lipid_species.txt", sep="\t")

# get reactions with lipid species -- lipid reactions

lipid_species_ID<-as.character(lipid_species$Metabolite)

lipid_species_rxn<-rxn_yeast[(rxn_yeast$Metabolite %in% lipid_species_ID),]
lipid_species_rxn<-lipid_species_rxn[!duplicated(lipid_species_rxn$ReactionID), ]

# get all species in lipid reactions

lipid_rxn_ID<-as.character(lipid_species_rxn$ReactionID)

lipid_rxn_StoiCoeff<-rxn_yeast[(rxn_yeast$ReactionID %in% lipid_rxn_ID),]
other_rxn_StoiCoeff<-rxn_yeast[!(rxn_yeast$ReactionID %in% lipid_rxn_ID),]

# build stoichiometric matrix

library("reshape2", lib.loc="D:/Dropbox/Eigene Dateien/Dokumente/R/win-library/3.0")

lipid<-lipid_rxn_StoiCoeff[!is.na(lipid_rxn_StoiCoeff$StoiCoef),]
S_lipid<-acast(lipid, Metabolite~ReactionID, value.var="StoiCoef")
S_lipid[is.na(S_lipid)]<-0

other<-other_rxn_StoiCoeff[!is.na(other_rxn_StoiCoeff$StoiCoef),]
S_other<-acast(other, Metabolite~ReactionID, value.var="StoiCoef")
S_other[is.na(S_other)]<-0

write.table(S_lipid, file="S_lipid.txt", sep="\t")
write.table(S_other, file="S_other.txt", sep="\t")

# reduce the stoichiometric matrix to reactions that carry flux without common species

S_lipid_reduced<-S_lipid[!(rownames(S_lipid) %in% rxn_yeast$Metabolite[grep('^H\\+|^H2O|^ATP |^ADP |^NAD|^phosphate|^ammonium', rxn_yeast$MetName)]),] # remove cofactors and common species 
S_lipid_reduced <- S_lipid_reduced[,-which(colSums(abs(S_lipid_reduced))==0)] # remove reactions that don't carry flux
S_other_reduced<-S_other[!(rownames(S_other) %in% rxn_yeast$Metabolite[grep('^H\\+|^H2O|^ATP |^ADP |^NAD|^phosphate|^ammonium', rxn_yeast$MetName)]),] # remove cofactors and common species 
S_other_reduced <- S_other_reduced[,-which(colSums(abs(S_other_reduced))==0)] # remove reactions that don't carry flux

# get common species in both subsets

common_species<-Reduce(intersect, list(row.names(S_lipid_reduced), row.names(S_other_reduced)))
common_species_rxn<-rxn_yeast[(rxn_yeast$Metabolite %in% common_species),]

# bipartite plot (igraph.incidence) of stoichiometric matrices in igraph

library(igraph)

lipid_reduced_graph <- graph.incidence(S_lipid_reduced, directed = T)
other_reduced_graph <- graph.incidence(S_other_reduced, directed = T)

# plot the bipartite plots

E(lipid_reduced_graph)$width = 0.1
E(lipid_reduced_graph)$weight = 0.1
E(lipid_reduced_graph)$color = "gray25"
V(lipid_reduced_graph)$label.cex = 1
V(lipid_reduced_graph)$label = clusters(lipid_reduced_graph)$membership
V(lipid_reduced_graph)$size = 2
V(lipid_reduced_graph)$color = ifelse(substr(V(lipid_reduced_graph)$name, 1, 1) == "s", "RED", "BLUE")

plot.igraph(lipid_reduced_graph, layout=layout.kamada.kawai, edge.arrow.size=0.1)
plot.igraph(lipid_reduced_graph, layout=layout.fruchterman.reingold, edge.arrow.size=0.1)

tkplot(lipid_reduced_graph, layout=layout.fruchterman.reingold, edge.arrow.size=0.1)

# check cluster composition

V(lipid_reduced_graph)[clusters(lipid_reduced_graph)$membership==3]

# get species and their corresponding cluster membership

vertices_clusters<-get.data.frame(lipid_reduced_graph, what="vertices")
vertices_clusters$cluster<-clusters(lipid_reduced_graph)$membership


# reduce lipid_Reduced_graph further to solely include headgroup ractions (without fatty acids)

list_headgroup<-c("ethanolamine", "serine", "inositol", "choline", "PC", "PE", "sphingo", "sphinganine", "sphingosine", "PS", "PI", "sterol", "TAG", "DAG", "diglyceride", "triglyceride", "PG", "CL", "phosphatidyl", "diacylglycerol", "triacylglycerol", "PA", "phosphatidate", "ceramide", "M(IP)2C", "2C")
lipid_headgroup<-lipid_rxn_StoiCoeff[which(grepl(paste(list_headgroup,collapse="|"), lipid_rxn_StoiCoeff$Reaction)),] # get headgroup reactions
lipid_headgroup<-lipid_headgroup[-which(grepl("isa", lipid_headgroup$Reaction)),] # remove "isa" pseudoreactions
lipid_headgroup<-lipid_headgroup[-which(grepl("exchange", lipid_headgroup$Compartment)),] # remove any exchange reactions
lipid_headgroup_reaction<-unique(lipid_headgroup$ReactionID)
lipid_headgroup_StoiCoeff<-rxn_yeast[(rxn_yeast$ReactionID %in% lipid_headgroup_reaction),]
lipid_headgroup_metabolites<-lipid_headgroup_StoiCoeff[!is.na(lipid_headgroup_StoiCoeff$StoiCoef),]
S_headgroup<-acast(lipid_headgroup_metabolites, Metabolite~ReactionID, value.var="StoiCoef")
S_headgroup[is.na(S_headgroup)]<-0
S_headgroup_reduced<-S_headgroup[!(rownames(S_headgroup) %in% rxn_yeast$Metabolite[grep('^H\\+|^H2O|^ATP |^ADP |^NAD|^phosphate|^ammonium', rxn_yeast$MetName)]),] # remove cofactors and common species
S_headgroup_reduced<-S_headgroup_reduced[,which(colSums(abs(S_headgroup_reduced))>0)]
headgroup_graph <- graph.incidence(S_headgroup_reduced, directed = T)

vertices_clusters_headgroup<-get.data.frame(headgroup_graph, what="vertices")
vertices_clusters_headgroup$cluster<-clusters(headgroup_graph)$membership
headgroup_speciesID<-vertices_clusters_headgroup[which(grepl("s_", vertices_clusters_headgroup$name)),]
headgroup_species<-lipid_headgroup_StoiCoeff[lipid_headgroup_StoiCoeff$Metabolite %in% headgroup_speciesID$name,]

headgroup_metabolites<-headgroup_species[!duplicated(headgroup_species$Metabolite), ]
headgroup_reactions<-headgroup_species[!duplicated(headgroup_species$ReactionID), ]

write.table(vertices_clusters_headgroup, "headgroup_clusters.csv", sep=",")
write.table(headgroup_metabolites, "headgroup_metabolites.csv", sep=",")
write.table(headgroup_reactions, "headgroup_reactions.csv", sep=",")

E(headgroup_graph)$width = 0.1
E(headgroup_graph)$weight = 0.1
E(headgroup_graph)$color = "gray25"
V(headgroup_graph)$label.cex = 1
V(headgroup_graph)$label = clusters(headgroup_graph)$membership
V(headgroup_graph)$size = 2
V(headgroup_graph)$color = ifelse(substr(V(headgroup_graph)$name, 1, 1) == "s", "RED", "BLUE")

tkplot(headgroup_graph, layout=layout.fruchterman.reingold, edge.arrow.size=0.1)

plot.igraph(headgroup_graph, layout=layout.kamada.kawai, edge.arrow.size=0.1)
plot.igraph(headgroup_graph, layout=layout.fruchterman.reingold, edge.arrow.size=0.1)
