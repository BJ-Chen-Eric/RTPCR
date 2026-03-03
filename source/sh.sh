Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260211jar319-2.csv,/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260211jar319.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-name output_20260211jar319_merge \
  --color-family blue \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 \
  --perform-calibration true \
  --calibration-sample M \
  --calibration-gene 5.8s \
  --calibration-time 0h \
  --calibration-rep 1 \
  --remove-calibration-sample true



Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260212jartub.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260211jarin tub.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260212jartub_merge' \
  --color-family blue \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 



Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260211jarin tub.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260212jartub.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260211jartub_merge' \
  --color-family blue \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 



Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260212cdk.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tcp43myb.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "MYB138,TCP44,TCP43" \
  --subset-topn-rules "s:MYB138:1:5,s:TCP44:1:5,s:TCP44:4:-5,l:TCP44:4:-5,l:TCP44:1:-5,l:TCP44:2:-5" \
  --subset-topn-rules-control "l:TCP44:0:4" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260216_target' \
  --color-family 'blue,green,pink' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 


Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260212cdk.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tcp43myb.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "CDC1,CDK3" \
  --subset-topn-rules "L:CDC1:0:-5" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260216_CDK' \
  --color-family 'blue,green,pink' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 


Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260214cdk_merge0205' \
  --color-family 'blue,green,pink' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 




Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251210hormone_dup1_biorep.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251211hormone_dup2_biorep.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251212hormone_dup3_biorep.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260212hormone_check' \
  --color-family 'blue,green,pink' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 




####### 0219 Pri319 
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260219premir319.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260209pri319_tub.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260219premir319' \
  --color-family 'blue,green,pink' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 




####### 0221 20260221arf5.csv, 20260221arf6 saur.csv
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf5.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf6 saur.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 1h \
  --target-genes "pin3,arf5,arf6,saur32" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260221arf5_arf6_saur' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 


####### 0221 strain compare (by_strain script)
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_strain_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221opr78sap_avg2rows.csv" \
  --reference-strain b73 \
  --subset-topn-rules-control "b73:mir319:0:-3" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260221opr78sap_by_strain' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 



####### 0221 20260221arf5.csv, 20260221arf6 saur.csv

Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf5.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf6 saur.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260223arf.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "pin3,arf5,arf6,saur32" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260223arf5_arf6_saur_merge' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 



Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf5.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260223arf.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "arf5,arf6" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/test' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 


####### 0225 volumn compare (by_strain script) 1
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_strain_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja.csv" \
  --reference-strain mock \
  # --subset-topn-rules-control "b73:mir319:0:-3" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260225ja_by_volumn' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 



####### 0225 volumn compare (by_strain script) 1+2
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_strain_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja-2.csv" \
  --reference-strain mock \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260225ja_by_volumn_1+2' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 


Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_strain_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja-2.csv" \
  --reference-strain mock \
  --subset-topn-rules "0.1um:mir319:0:-4,0.5um:mir319:0:4,1um:mir319:0:4,10um:mir319:0:4" \
  --sample-order "0.1um,0.5um,1um,10um" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260225ja_by_volumn_1+2' \
  --color-family 'orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 




Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_strain_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja_avg_ct_by_row.csv" \
  --reference-strain mock \
  --subset-topn-rules "0.1um:mir319:0:-3" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260225ja_avg' \
  --color-family 'green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 




####### 0228 20260228sorganmir319
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260228sorganmir319_2.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260228sorganmir319.csv" \
  --reference-strain mock \
  --subset-topn-rules "s:mir319:1:-5" \
  --subset-topn-rules-control "l:mir319:0:-4" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260228sorganmir319' \
  --color-family 'teal' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 


##### figure update
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified_figure_update.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf5.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260223arf.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "arf5,arf6" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/test' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 




Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251219localmir319.csv,,/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251223mir159lical2.csv" \
  --reference-strain mock \
  --subset-topn-rules "B73:mir319:1:-5" \
  --subset-topn-rules-control "B73:mir319:0:5" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260303_1219_test' \
  --color-family 'purple' \
  --color-seed 929 \
  --color-count 12 \
  --color-index 1 

/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251219localmir319.csv,
,/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251223mir159lical2.csv