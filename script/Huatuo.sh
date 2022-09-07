#---------- define user input parameter -----------#

run_script_folder=$(pwd)
input_variant_file=${run_script_folder}"/example_data/input_variant.txt"
anno_file=${run_script_folder}"/data/HCL_Cell_Cluster_Annatation.txt"
XGBoost_modelslist=${run_script_folder}"/data/HCL_XGBoost_modelslist"

while [ $# -ge 2 ]; do
        case "$1" in
                --input_variant_file)
                        input_variant_file=$2; shift 2;;
                --anno_file)
                        anno_file=$2; shift 2;;
                --XGBoost_modelslist)
                        XGBoost_modelslist=$2; shift 2;;
                *)
                        echo "unknown input parameter!"; exit 1;;
        esac
done

#------------ define internal path --------------#

result_directory=${run_script_folder}"/result/"
CNN_result_directory=${run_script_folder}"/result/de_novo_predictions/"
Landscape_ieQTL_file=${run_script_folder}"/data/Landscape_ieQTL.txt"
Landscape_ieQTL_Linked_SNPs_file=${run_script_folder}"/data/Landscape_ieQTL_Linked_SNPs.RData"
coordinate_file=${run_script_folder}"/resource/hg19ToHg38.over.chain.gz"
hg19_fa=${run_script_folder}"/resource/hg19.fa"
h19_geneanno_file=${run_script_folder}"/resource/geneanno.csv"
beluga_model=${run_script_folder}"/resource/deepsea.beluga.pth"

mkdir $result_directory
mkdir $CNN_result_directory

#------------- variant predict -----------#

echo 'Computing the chromatin effects of the variants...'
Rscript ${run_script_folder}"/script/Prepare_predicted_vcf.r" --input_variant_file $input_variant_file\
   --h19_geneanno_file $h19_geneanno_file \
   --CNN_result_directory $CNN_result_directory

cd ${CNN_result_directory}
python ${run_script_folder}"/script/chromatin.py" predicted.vcf \
       --hg19_fa $hg19_fa \
       --deepsea $beluga_model 

echo 'Computing variant effects for each cell type...'
cd ${run_script_folder}
python ${run_script_folder}"/script/predict.py" --coorFile ${CNN_result_directory}predicted.vcf \
            --snpEffectFilePattern ${CNN_result_directory}predicted.vcf.shift_SHIFT.diff.h5 \
            --modelList $XGBoost_modelslist --output ${CNN_result_directory}prediction_output_no_gene.csv

#--------------- inferring cell-type-specific gene regulation ------------#

echo 'Inferring cell-type-specific gene regulation...'

Rscript ${run_script_folder}"/script/Inferring_gene_regulation.r" --CNN_result_directory $CNN_result_directory \
  --result_directory $result_directory \
  --Landscape_ieQTL_file $Landscape_ieQTL_file \
  --coordinate_file $coordinate_file \
  --Landscape_ieQTL_Linked_SNPs_file $Landscape_ieQTL_Linked_SNPs_file \
  --anno_file $anno_file









