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
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260211jarin tub.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/20260211jartub' \
  --color-family blue \
  --color-seed 123 \
  --color-count 12 \
  --color-index 1 

