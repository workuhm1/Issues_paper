preprocess_pdf <- function(pdf_txt_raw = NULL) {
  
  pdf_page2 <- data.frame(a=numeric(), b=character())
  count <- 0
  
  for (txtin in pdf_txt_raw) {
    
    count <- count + 1
    sys_name <- tolower(Sys.info()["sysname"])
    
    # if (!(str_detect(def_OS, "windows"))) {
    if (!stringr::str_detect(sys_name, "windows")) {
      txtin_splt <- unlist(strsplit(as.character(txtin), "\n"))
    }
    else {
      txtin_splt <- unlist(strsplit(as.character(txtin), "\r\n"))
    }
    
    txtin_splt <- cbind(count, txtin_splt)
    
    pdf_page2 <- rbind(pdf_page2, txtin_splt)
  }
  
  pdf_page2 <- as.data.frame(pdf_page2)
  colnames(pdf_page2) <- c("page_nbr", "pdf_txt")
  
  return(pdf_page2)
  
}