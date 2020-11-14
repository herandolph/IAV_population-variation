#!/bin/bash

module load java/1.8

CONDITIONS=(flu NI)
PSEUDOBULK=(pseudosums_QN)
PCS=(3 6)
TEMP_DIRS=(3expPC_Flu 6expPC_NI)
OUT_DIRS=(filtered_3expPCs_propBimodal filtered_6expPCs_propBimodal)
META=(individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt)
COUNTS=(corrected_expression.txt)
GENE_POS=(GRCh38.92_gene_positions.txt)
CELLTYPE=(B)

LEN=${#CONDITIONS[@]}

for (( NUM=0; NUM<$LEN; NUM++ ))
    do
    	
       CONDITION=${CONDITIONS[$NUM]}
       PC=${PCS[$NUM]}
       TEMP_DIR=${TEMP_DIRS[$NUM]}
       OUT_DIR=${OUT_DIRS[$NUM]}
       

       echo " ************************************** "
       echo " BEGIN MATRIXEQTL PIPELINE FOR: $PSEUDOBULK "
       echo " CONDITION = $CONDITION;"
       echo " EXP PC TO REGRESS = $PC;"
       echo " TEMP DIRECTORY = $TEMP_DIR;"
       echo " OUTPUT DIRECTORY = $OUT_DIR;"
       echo " ************************************** "
       sbatch eQTL_job_specs.sbatch $CONDITION $PSEUDOBULK $PC $TEMP_DIR $OUT_DIR $META $COUNTS $GENE_POS $CELLTYPE &
    done
wait

echo "ALL JOBS SUMBITTED AT: `date`"

## EOF ##
