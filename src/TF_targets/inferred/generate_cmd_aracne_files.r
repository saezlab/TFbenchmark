rm(list = ls())
home = '/Volumes/GoogleDrive/My\\ Drive/projects/TFbenchmark/code/TF_targets/inferredGTEx/'
setwd(home)


# Load expression dataset
exp = get(load('~/data/GTEx_expressionatlas/raw/GTex_v6_atlas.genes.voom.batchcor.merged.rdata'))


# Load samples annotation
load('~/data/GTEx_expressionatlas/annot/E-MTAB-5214.sdrf.rdata')


# Define tissues
tissues = sample_annotation$Comment.histological.type.[ sample_annotation$Source.Name %in% colnames(exp) ]
tissues = sort(unique(tissues))


# Define command lines
step_00 = "echo \"--------------------\"MyTissueDirectory\"--------------------\"; echo; echo"
step_0 = "mkdir -p ~/tmp/inferredGTExMyTissueDirectory/"
step_1 = "java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/MyTissueDirectory.txt -o ~/tmp/inferredGTExMyTissueDirectory/ --tfs /Volumes/GoogleDrive/My\\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold"
step_2 = "for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/MyTissueDirectory.txt -o ~/tmp/inferredGTExMyTissueDirectory/ --tfs /Volumes/GoogleDrive/My\\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done"
step_3 = "java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExMyTissueDirectory/ --consolidate"


CMD = list()
for (ti in tissues) {
  ti = gsub(' ', '_', ti)
  CMD = append(CMD, gsub('MyTissueDirectory', ti, step_00))
  CMD = append(CMD, gsub('MyTissueDirectory', ti, step_0))
  CMD = append(CMD, gsub('MyTissueDirectory', ti, step_1))
  CMD = append(CMD, gsub('MyTissueDirectory', ti, step_2))
  CMD = append(x = CMD, values = gsub('MyTissueDirectory', ti, step_3))
}
write.table(unlist(CMD), file = 'run_aracne.cmd', quote = F, row.names = F, col.names = F)
