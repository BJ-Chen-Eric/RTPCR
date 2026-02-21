Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
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



Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260212jartub.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260211jarin tub.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260212jartub_merge' \
  --color-family blue \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 



Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260211jarin tub.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260212jartub.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260211jartub_merge' \
  --color-family blue \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 



Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
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


Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
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


Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260214cdk_merge0205' \
  --color-family 'blue,green,pink' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 




Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
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
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
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
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_pct_strain_source_unified.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf5.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf6 saur.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "pin3,arf5,arf6,saur32" \
  --subset-topn-rules "S:arf6:4:5" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260221arf5_arf6_saur' \
  --color-family 'blue,green,pink,orange' \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 
