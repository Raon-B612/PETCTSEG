##-- 00. initialize --##
rm(list=ls())
source('~/CONFIG_EMMA.R')
source('~/PETCTSEG/scripts/CONFIG_PETCTSEG.R')

# set working directory
src_dir <- path('C:/Temp/PETCTSRC/')
mkdir(src_dir)

# quick update segmentation database
seg_db <- read_csv(path(share_dir, 'seg_db.csv'), show_col_types = FALSE)
seg_db[1,] |> select(pdate, qdate)

#-- 00. preprocessing --#
dcm_df <- read_dicom_dir('F:/PETCTDCM/')
dcm2nii_ct(dcm_df)
dcm2nii_pet(dcm_df)

#-- 01. resample PET to CT space --#
rename_lx(src_dir)
move_deprecated(src_dir)
resample_suv_to_ct(src_dir)

# -- 02. segmentation --#
prepare_segmentation(src_dir, task = 'all')
ts <- read_csv(path(prep_dir, 'tseg_prep.csv'),  show_col_types = FALSE)
ms <- read_csv(path(prep_dir, 'moose_prep.csv'), show_col_types = FALSE)
# segment_ct(src_dir)

#-- 03. roi_measurement and summary --#
measure_roi(src_dir, overwrite = TRUE)
summarize_subject_measures(overwrite = TRUE)

#-- 04. distribute_processed_files --#
distribute_processed_files(gz_dir)

#-- 05. update segmentation database --#
update_seg_db()
