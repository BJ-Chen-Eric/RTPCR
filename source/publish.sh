##### figure update demo
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified_figure_update.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf5.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260223arf.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_remake/test' 
  


####### 0206 L/S mir319 
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified_figure_update.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260206mir319b73.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260206mir3192.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --subset-topn-rules "s:mir319:1:-3" \
  --subset-topn-rules-control "s:mir319:0:4" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_publish/0206mir319_LS' 
  


####### 0216 L/S CDK3
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified_figure_update.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260212cdk.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tcp43myb.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "CDK3" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_publish/20260216_CDK' 
  



####### 0221 L/S 20260221arf5.csv, 20260221arf6 saur.csv
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified_figure_update.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf5.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260214myb138.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260205tub_tcp44.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260221arf6 saur.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260223arf.csv" \
  --analysis-mode multi \
  --multi-compare-style all_time \
  --control-time 0h \
  --target-genes "arf5,arf6" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_publish/20260221arf' 





####### 0221 L/S 20260221arf5.csv, 20260221arf6 saur.csv
# Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_strain_unified.R \
#   --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251210hormone_dup1_biorep.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251211hormone_dup2_biorep.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20251212hormone_dup3_biorep.csv" \
#    --reference-strain c \
#   --analysis-mode multi \
#   --multi-compare-style all_time \
#   --control-time 0h \
#   --target-genes "arf5,arf6" \
#   --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_publish/20260221arf' 


####### JA mir319
# /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline/rt_pct_strain_source_1227_hor_0227.R



####### 0228 20260228sorganmir319
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_time_unified_figure_update.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260228sorganmir319_2.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260228sorganmir319.csv" \
  --reference-strain mock \
  --subset-topn-rules "s:mir319:1:-5" \
  --subset-topn-rules-control "l:mir319:0:-4" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_publish/20260228sorganmir319'
  

####### 0225 volumn compare (by_strain script) 1+2
Rscript /Users/minieric/Desktop/projects/Jean_maize/Jean/pipeline_v2/work_pp/rt_strain_unified_figure_update.R \
  --files "/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja.csv, /Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/raw/20260225ja-2.csv" \
  --reference-strain mock \
  --subset-topn-rules "0.1um:mir319:0:-4,0.5um:mir319:0:4,1um:mir319:0:4,10um:mir319:0:4" \
  --sample-order "0.1um,0.5um,1um,10um" \
  --out-dir '/Users/minieric/Desktop/projects/Jean_maize/Jean/RTpcr/out_publish/20260225ja_by_volumn_1+2' 