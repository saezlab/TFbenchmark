rm(list = ls())
home = '/Volumes/GoogleDrive/My Drive/projects/TFbenchmark/'
setwd(home)
source('code/utils.r')



TFrole_genesets = list()


# TF mode_of_regulation 
activators = load_activators()
dualActRep = load_dualActRep()
repressors = load_repressors()
TFrole_genesets[['mode_of_regulation']][['activators']] = activators
TFrole_genesets[['mode_of_regulation']][['repressors']] = repressors
TFrole_genesets[['mode_of_regulation']][['dual']] = dualActRep



# TF DNA_binding_mode
complexes = load_complexes()
TFrole_genesets[['DNA_binding_mode']][['Obligate_heteromer']] = unique(subset(complexes, complex_type == 'heterodimer')$gene_symbol)
TFrole_genesets[['DNA_binding_mode']][['Monomer_or_homomultimer']] = unique(subset(complexes, complex_type != 'heterodimer')$gene_symbol)
# from lambet2018
complexes = read.csv('data/TF_info/regulation_type/lambert_cell2018/TableS1.csv', stringsAsFactors = F)
TFrole_genesets[['DNA_binding_mode']][['Obligate_heteromer']] = unique(c(TFrole_genesets[['DNA_binding_mode']][['Obligate_heteromer']], subset(complexes, X.2 == '2 Obligate heteromer')$X))
TFrole_genesets[['DNA_binding_mode']][['Monomer_or_homomultimer']] = unique(c(TFrole_genesets[['DNA_binding_mode']][['Monomer_or_homomultimer']] , subset(complexes, X.2 == '1 Monomer or homomultimer')$X))
contradictory = intersect(TFrole_genesets[['DNA_binding_mode']][['Monomer_or_homomultimer']], 
                        TFrole_genesets[['DNA_binding_mode']][['Obligate_heteromer']])
TFrole_genesets[['DNA_binding_mode']][['Obligate_heteromer']] = setdiff(TFrole_genesets[['DNA_binding_mode']][['Obligate_heteromer']],
                                                                             contradictory)
TFrole_genesets[['DNA_binding_mode']][['Monomer_or_homomultimer']] = setdiff(TFrole_genesets[['DNA_binding_mode']][['Monomer_or_homomultimer']],
                                                                             contradictory)



# DNA_binding_specificity from lambet2018
TFrole_genesets[['DNA_binding_specificity']][['Low_specificity_DNA-binding_protein']] = unique(subset(complexes, X.2 == '3 Low specificity DNA-binding protein')$X)
TFrole_genesets[['DNA_binding_specificity']][['Not_a_DNA_binding_protein']] = unique(subset(complexes, X.2 == '4 Not a DNA binding protein')$X)
TFrole_genesets[['DNA_binding_specificity']][['High_specificity_DNA-binding_protein']] = unique(subset(complexes, X.2 %in% c('1 Monomer or homomultimer', '2 Obligate heteromer'))$X)





# 
# # TF chromatin regulation
CHrole = load_CHregulation()
TFrole_genesets[['chromatin_regulation_mode']][['settler']] = unique(subset(CHrole, Chromatin_Opening_Type == 'SETTLER')$Approved_symbol)
TFrole_genesets[['chromatin_regulation_mode']][['migrant']] = unique(subset(CHrole, Chromatin_Opening_Type %in% c('NEGATIVE_MIGRANT', 'POSITIVE_MIGRANT'))$Approved_symbol)
TFrole_genesets[['chromatin_regulation_mode']][['pioneer']] = unique(subset(CHrole, Chromatin_Opening_Type == 'PIONEER')$Approved_symbol)




# TF_class
TFclassification = load_TFclass_classes()
for (gr in unique(TFclassification$superclass_name) )
  TFrole_genesets[['TF_superclass']][[gr]] = unique(subset(TFclassification, superclass_name == gr)$name)
for (gr in unique(TFclassification$class_name) )
  TFrole_genesets[['TF_class']][[gr]] = unique(subset(TFclassification, class_name == gr)$name)
# for (gr in unique(TFclassification$family_name) )
#   TFrole_genesets[['TF_class']][[gr]] = unique(subset(TFclassification, family_name == gr)$name)




# # TF tissue_of_expression
TF_x_tissue = load_TFtissues()
TF_x_tissue_freq = table(TF_x_tissue$V1) / max(table(TF_x_tissue$V1))
TFrole_genesets[['tissue_of_expression']][['tissue-specific']] = names(which(TF_x_tissue_freq<0.1))
TFrole_genesets[['tissue_of_expression']][['no_tissue-specific']] = names(which(TF_x_tissue_freq>0.9))




# Save
save(TFrole_genesets, file = 'data/TF_info/TFrole_genesets.rdata')
