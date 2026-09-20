##-- 00. initialize --##
rm(list=ls())
source('~/PETCTSEG/config_petctseg.R')

df_db <- dir_ls(gz_all, type = 'directory', recurse = FALSE,regexp = '_\\d{8}' ) |>
  tibble(db_dir = _) |>
  mutate(
    subjid_long = basename(db_dir),
    subjid = str_remove(subjid_long, '_\\d{8}'),
    status = 'segmented',
    qdate = Sys.Date() |> ymd()
  ) |>
  arrange(subjid)

df_summary <- dir_ls(summary_dir, type = 'file', recurse = TRUE, regexp = '\\d{10,13}\\.csv$') |>
  tibble(summary = _) |>
  mutate(
    subjid = str_remove_all(basename(summary), 'summary_|\\.csv'),
    status = 'summarized'
    )

df_measure <- dir_ls(measure_dir, type = 'directory', recurse = TRUE, regexp = '_\\d{8}') |>
  tibble(csv_dir = _) |>
  mutate(subjid_long = basename(csv_dir),
         subjid = str_remove(subjid_long, '_\\d{8}'),
         accn = str_extract(subjid_long, '\\d{10,13}'),
         status = 'measured'
  )

df_proc <- dir_ls(gz_dir, type = 'directory', recurse = TRUE, regexp = '_\\d{8}') |>
  tibble(gz_dir = _) |>
  mutate(
    subjid_long = basename(gz_dir),
    subjid = str_remove(subjid_long, '_\\d{8}'),
    status = 'processing'
    )

df_arx <- dir_ls('Z:/PETCTSEG/', type = 'directory', recurse = FALSE, regexp = '[BMPVX]$') |>
  dir_ls(type = 'directory', recurse = FALSE, regexp = '_\\d{8}') |>
  tibble(arx_dir = _) |>
  mutate(
    subjid_long = basename(arx_dir),
    subjid = str_remove(subjid_long, '_\\d{8}')
  )

df_all <- bind_rows(df_db, df_proc, df_measure, df_summary) |>
  mutate(qdate = Sys.Date() |> ymd())

write_csv(df_db,  path(share_dir, 'seg_db_short.csv'))
write_csv(df_all, path(share_dir, 'df_all.csv'))

stop()
# rm_dl <- filter(df_all, subjid == 'ENM1211008454') |>
#   select(db_dir, csv_dir, gz_dir, summary) |>
#   unlist() |>
#   unname() %>%
#   keep(!is.na(.))
# 
# print(rm_dl)
# file_delete(rm_dl)

stop()
# rm_dl <- left_join(df_proc, df_db, by = 'subjid_long') |>
#   filter(!is.na(db_dir))
# rm_dl
df_proc
stop()
rm_dl <- NULL
for (i in 1:nrow(df_proc)) {
  # i <- 1
  x <- df_proc[i,]
  files <- dir_ls(x$gz_dir, type = 'file', recurse = FALSE)
  zero_files <- files[file_size(files) == 0]
  if (length(zero_files != 0)) {
    rm_dl <- rbind(rm_dl, x)
    }
}
rm_dl
write.csv(rm_dl, path(share_dir, 'zero_files.csv'))

stop()
rm_dl <- NULL
rm_dl <- left_join(df_arx, df_db, by = 'subjid') |>
  filter(is.na(qdate))
rm_dl
# rm_dl |> pull(arx_dir) |> file_delete()

stop()
left_join(df_proc, seg_db, by='subjid_long') |>
  filter(is.na(db_dir))
seg_db
