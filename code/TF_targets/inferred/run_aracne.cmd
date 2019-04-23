echo "--------------------"Adipose_Tissue"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExAdipose_Tissue/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Adipose_Tissue.txt -o ~/tmp/inferredGTExAdipose_Tissue/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Adipose_Tissue.txt -o ~/tmp/inferredGTExAdipose_Tissue/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExAdipose_Tissue/ --consolidate
echo "--------------------"Adrenal_Gland"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExAdrenal_Gland/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Adrenal_Gland.txt -o ~/tmp/inferredGTExAdrenal_Gland/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Adrenal_Gland.txt -o ~/tmp/inferredGTExAdrenal_Gland/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExAdrenal_Gland/ --consolidate
echo "--------------------"Bladder"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExBladder/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Bladder.txt -o ~/tmp/inferredGTExBladder/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Bladder.txt -o ~/tmp/inferredGTExBladder/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExBladder/ --consolidate
echo "--------------------"Blood"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExBlood/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Blood.txt -o ~/tmp/inferredGTExBlood/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Blood.txt -o ~/tmp/inferredGTExBlood/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExBlood/ --consolidate
echo "--------------------"Blood_Vessel"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExBlood_Vessel/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Blood_Vessel.txt -o ~/tmp/inferredGTExBlood_Vessel/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Blood_Vessel.txt -o ~/tmp/inferredGTExBlood_Vessel/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExBlood_Vessel/ --consolidate
echo "--------------------"Brain"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExBrain/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Brain.txt -o ~/tmp/inferredGTExBrain/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Brain.txt -o ~/tmp/inferredGTExBrain/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExBrain/ --consolidate
echo "--------------------"Breast"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExBreast/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Breast.txt -o ~/tmp/inferredGTExBreast/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Breast.txt -o ~/tmp/inferredGTExBreast/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExBreast/ --consolidate
echo "--------------------"Cervix_Uteri"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExCervix_Uteri/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Cervix_Uteri.txt -o ~/tmp/inferredGTExCervix_Uteri/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Cervix_Uteri.txt -o ~/tmp/inferredGTExCervix_Uteri/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExCervix_Uteri/ --consolidate
echo "--------------------"Colon"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExColon/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Colon.txt -o ~/tmp/inferredGTExColon/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Colon.txt -o ~/tmp/inferredGTExColon/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExColon/ --consolidate
echo "--------------------"Esophagus"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExEsophagus/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Esophagus.txt -o ~/tmp/inferredGTExEsophagus/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Esophagus.txt -o ~/tmp/inferredGTExEsophagus/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExEsophagus/ --consolidate
echo "--------------------"Fallopian_Tube"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExFallopian_Tube/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Fallopian_Tube.txt -o ~/tmp/inferredGTExFallopian_Tube/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Fallopian_Tube.txt -o ~/tmp/inferredGTExFallopian_Tube/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExFallopian_Tube/ --consolidate
echo "--------------------"Heart"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExHeart/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Heart.txt -o ~/tmp/inferredGTExHeart/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Heart.txt -o ~/tmp/inferredGTExHeart/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExHeart/ --consolidate
echo "--------------------"Kidney"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExKidney/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Kidney.txt -o ~/tmp/inferredGTExKidney/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Kidney.txt -o ~/tmp/inferredGTExKidney/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExKidney/ --consolidate
echo "--------------------"Liver"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExLiver/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Liver.txt -o ~/tmp/inferredGTExLiver/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Liver.txt -o ~/tmp/inferredGTExLiver/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExLiver/ --consolidate
echo "--------------------"Lung"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExLung/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Lung.txt -o ~/tmp/inferredGTExLung/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Lung.txt -o ~/tmp/inferredGTExLung/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExLung/ --consolidate
echo "--------------------"Muscle"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExMuscle/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Muscle.txt -o ~/tmp/inferredGTExMuscle/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Muscle.txt -o ~/tmp/inferredGTExMuscle/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExMuscle/ --consolidate
echo "--------------------"Nerve"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExNerve/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Nerve.txt -o ~/tmp/inferredGTExNerve/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Nerve.txt -o ~/tmp/inferredGTExNerve/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExNerve/ --consolidate
echo "--------------------"Ovary"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExOvary/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Ovary.txt -o ~/tmp/inferredGTExOvary/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Ovary.txt -o ~/tmp/inferredGTExOvary/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExOvary/ --consolidate
echo "--------------------"Pancreas"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExPancreas/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Pancreas.txt -o ~/tmp/inferredGTExPancreas/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Pancreas.txt -o ~/tmp/inferredGTExPancreas/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExPancreas/ --consolidate
echo "--------------------"Pituitary"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExPituitary/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Pituitary.txt -o ~/tmp/inferredGTExPituitary/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Pituitary.txt -o ~/tmp/inferredGTExPituitary/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExPituitary/ --consolidate
echo "--------------------"Prostate"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExProstate/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Prostate.txt -o ~/tmp/inferredGTExProstate/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Prostate.txt -o ~/tmp/inferredGTExProstate/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExProstate/ --consolidate
echo "--------------------"Salivary_Gland"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExSalivary_Gland/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Salivary_Gland.txt -o ~/tmp/inferredGTExSalivary_Gland/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Salivary_Gland.txt -o ~/tmp/inferredGTExSalivary_Gland/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExSalivary_Gland/ --consolidate
echo "--------------------"Skin"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExSkin/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Skin.txt -o ~/tmp/inferredGTExSkin/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Skin.txt -o ~/tmp/inferredGTExSkin/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExSkin/ --consolidate
echo "--------------------"Small_Intestine"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExSmall_Intestine/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Small_Intestine.txt -o ~/tmp/inferredGTExSmall_Intestine/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Small_Intestine.txt -o ~/tmp/inferredGTExSmall_Intestine/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExSmall_Intestine/ --consolidate
echo "--------------------"Spleen"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExSpleen/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Spleen.txt -o ~/tmp/inferredGTExSpleen/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Spleen.txt -o ~/tmp/inferredGTExSpleen/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExSpleen/ --consolidate
echo "--------------------"Stomach"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExStomach/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Stomach.txt -o ~/tmp/inferredGTExStomach/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Stomach.txt -o ~/tmp/inferredGTExStomach/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExStomach/ --consolidate
echo "--------------------"Testis"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExTestis/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Testis.txt -o ~/tmp/inferredGTExTestis/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Testis.txt -o ~/tmp/inferredGTExTestis/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExTestis/ --consolidate
echo "--------------------"Thyroid"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExThyroid/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Thyroid.txt -o ~/tmp/inferredGTExThyroid/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Thyroid.txt -o ~/tmp/inferredGTExThyroid/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExThyroid/ --consolidate
echo "--------------------"Uterus"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExUterus/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Uterus.txt -o ~/tmp/inferredGTExUterus/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Uterus.txt -o ~/tmp/inferredGTExUterus/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExUterus/ --consolidate
echo "--------------------"Vagina"--------------------"; echo; echo
mkdir -p ~/tmp/inferredGTExVagina/
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Vagina.txt -o ~/tmp/inferredGTExVagina/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed 1 --calculateThreshold
for i in `seq 1 100`; do java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -e /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_target_sources/inferredGTEx/GTEx_tissue_matrix/Vagina.txt -o ~/tmp/inferredGTExVagina/ --tfs /Volumes/GoogleDrive/My\ Drive/projects/TFbenchmark/data/TF_census/TFClass/huTF_census.txt  --pvalue 1E-8 --seed $i; done
java -Xmx5G -jar /Users/luzgaral/software/ARACNe/Aracne.jar -o ~/tmp/inferredGTExVagina/ --consolidate
