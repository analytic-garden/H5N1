#' fasta2dataframe - read a GISAID EpiFlu fasta file to a dataframe
#'
#' @param fasta_file - a string containg that fasta file name
#' @param remove_ref - boolean. Should reference sequence be removed from dataframe
#' @param ref_row - the row number of teh reference sequence, default = 1. Ignored if remove_ref = FALSE
#' @param min_seq_length - sequences with few elements are removed
#'
#' @returns - a dataframe with 16 columns
#' "Influenza_type"   "Host"             "Country"          "ID"               "Year"             "Isolate_name"     "Isolate_ID"      
#' "Type"             "Clade"            "Collection_date"  "Segment"          "Lineage"          "DNA_Accession_no" "seq.name"        
#' "seq.text"         "seq_length"
#'
fasta2dataframe <- function(fasta_file, 
                            remove_ref = FALSE, 
                            ref_row = 1,
                            min_seq_length = 1700) {
  require(phylotools)  # for reading the FASTA file
  require(tidyverse)
  
  # a list of US states for converting atates to USA
  states <- c("Alabama", "Alaska", "Arizona",  "Arkansas", "California", "Colorado",       
              "Connecticut", "Delaware", "Florida","Georgia", "Hawaii",
              "Idaho","Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana",
              "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",      
              "Mississippi", "Missouri", "Montana",  "Nebraska", "Nevada",  
              "New_Hampshire", "New_Jersey", "New_Mexico", "New_York", "North_Carolina", "North_Dakota",
              "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode_Island", 
              "South_Carolina", "South_Dakota", "Tennessee", "Texas", "Utah", 
              "Vermont", "Virginia", "Washington", "West,Virginia", "Wisconsin", "Wyoming")
  
  df <- read.fasta(fasta_file)
  if(remove_ref) {
    df <- df %>% slice(-ref_row)
  }
  
  # split sequence  IDs and get rid of any sequence that doesn't have full information
  df <- df %>% separate_wider_delim(seq.name, 
                                    delim = "|", 
                                    names = c("Isolate_name", "Isolate_ID", "Type", "Clade", 
                                              "Collection_date", "Segment", "Lineage", "DNA_Accession_no"),
                                    too_many = "debug", too_few = "debug")
  df <- df %>% filter(`seq.name_ok`)
  
  df <- df %>% 
    separate_wider_delim(`Isolate_name`, 
                         delim = "/", 
                         names = c("Influenza_type", "Host", "Country", "ID", "Year"), 
                         too_many = "debug", too_few = "debug")
  df <- df %>% filter(`Isolate_name_ok`)
  
  # map states to USA and map the German entries to Germany. Foramt a few others.
  df <- df %>% mutate(Country = if_else(Country %in% states, "USA", Country))
  df <- df %>% mutate(Country = if_else(grepl("Germany", Country), "Germany", Country))
  df <- df %>% mutate(Country = if_else(grepl("Nordrhein-Westfalen", Country), "Germany", Country))
  df <- df %>% mutate(Country = if_else(Country %in% c("Sao_Paulo", "Rio_de_Janeiro"), "Brazil", Country))
  df <- df %>% mutate(Country = if_else(Country %in% c("Kagoshima", "Ishikawa", "Iwate", "Hokkaido", "Hiroshima"), "Japan", Country))
  
  # Combine common hosts
  df <- df %>% mutate(Host = str_to_upper(Host))
  df <- df %>% mutate(Host = if_else(grepl("DUCK", Host), "DUCK", Host))
  df <- df %>% mutate(Host = if_else(grepl("MALLARD", Host), "DUCK", Host))
  df <- df %>% mutate(Host = if_else(grepl("HEN", Host), "CHICKEN", Host))
  df <- df %>% mutate(Host = if_else(grepl("CAT", Host), "CAT", Host))
  df <- df %>% mutate(Host = if_else(grepl("FELINE", Host), "CAT", Host))
  df <- df %>% mutate(Host = if_else(grepl("GOOSE", Host), "GOOSE", Host))
  df <- df %>% mutate(Host = if_else(grepl("CROW", Host), "CROW", Host))
  df <- df %>% mutate(Host = if_else(grepl("FOX", Host), "FOX", Host))
  
  # filter sequence columns
  df <- df %>% filter(str_length(seq.text) >= min_seq_length) 
  df <- df %>% mutate(seq.text = str_to_upper(seq.text))
  df <- df %>% mutate(seq_length = str_length(seq.text))
  
  # Since we removed the rows with missing values, remove the debug columns
  df <- df %>% select(-c("Isolate_name_ok", "Isolate_name_pieces", "Isolate_name_remainder"))
  df <- df %>% select(-c("seq.name_ok", "seq.name_pieces", "seq.name_remainder"))
  
  return(df)
}