setwd(here::here())

read_metaG_metadata <- function(){

  import_data <- read_tsv("data/import_log.tsv",col_names = FALSE) %>% 
    dplyr::select(sample_id = "X1",
                  import_id = "X8") %>% 
    distinct()
  
  metadata <- read_tsv("data/sample_data/metadata.txt") %>% 
    left_join(import_data) %>% 
    filter(!is.na(import_id)) %>% 
    mutate(transplant_type = if_else(str_detect(transplant_type, "Na.ve"), "Naive", transplant_type),
           days_after_transplant = as.numeric(if_else(days_after_transplant == "N/A", "NA", days_after_transplant))) %>% 
    mutate(sample = import_id,
           gut_section = factor(gut_section, levels = c("Terminal ileum", "Cecum", "Transverse Colon", "Descending Colon"), ordered = TRUE),
           gut_section_simp = factor(gut_section, levels = c("Terminal ileum", "Cecum", "Transverse Colon", "Descending Colon"), ordered = TRUE),
           replicate = case_when(str_detect(researcher_sample_id, "A1") ~ "A",
                                 str_detect(researcher_sample_id, "A2") ~ "B",
                                 str_detect(researcher_sample_id, "A3") ~ "C",
                                 str_detect(researcher_sample_id, "S1") ~ "A",
                                 str_detect(researcher_sample_id, "S2") ~ "B",
                                 str_detect(researcher_sample_id, "S3") ~ "C",
                                 .default = "other"),
           days_after_transplant = case_when(replicate == "B" & days_after_transplant == 21 ~ "7",
                                             replicate == "B" & days_after_transplant == 7 ~ "21",
                                             transplant_type == "Naive" ~ "other",
                                             .default = as.character(days_after_transplant)),
           days_after_transplant = factor(days_after_transplant,levels = c("7","21","other"), ordered = TRUE),
           plot_categories = glue::glue("{transplant_type}; day {days_after_transplant}"),
           plot_categories = if_else(transplant_type == "Naive", "Naive", plot_categories),
           plot_categories = factor(plot_categories,
                                    levels = c("Naive","Syngeneic; day 7", "Syngeneic; day 21", "Allogeneic; day 7", "Allogeneic; day 21"), ordered = TRUE))

}