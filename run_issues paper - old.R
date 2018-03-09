
## load library
source("./load library.R")

## load supportive functions
source("./preprocess_pdf.R")

## Step 1. Get simulation results (from SEM article)
nonconvg_pdf <- pdftools::pdf_text(pdf = "./rawdata - pdf tables/Table 3_Nonconvergence.pdf")

## Step 2. Preprocess PDF file and get analysis dataset
## -- 2.1. Preprocess PDF file
nonconvg_pdf_prepro <- preprocess_pdf(pdf_txt_raw = nonconvg_pdf)

## add table row number
nonconvg_pdf_prepro2 <- nonconvg_pdf_prepro %>% 
  dplyr::filter(!row_number() %in% c(1:6, 17)) %>% 
  dplyr::mutate_if(is.factor, as.character) ## Convert data.frame columns from factors to characters (source: https://stackoverflow.com/questions/2851015/convert-data-frame-columns-from-factors-to-characters) 

## -- 2.2. Create analysis dataset
nonconvg_txt <- nonconvg_pdf_prepro2 %>% 
  dplyr::select(-page_nbr) %>% 
  dplyr::pull(pdf_txt) %>% 
  stringr::str_trim()
   
## Extract numbers
nonconvg_txt2 <- nonconvg_txt %>% 
  stringr::str_match_all(pattern = "[0-9]+") %>% 
  unlist() %>% 
  stringr::str_trim()


