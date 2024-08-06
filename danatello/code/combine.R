library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(optparse)

# Read table and add Genome name from the filename
read_tsv_name = function(filepath){
  filename = str_remove(basename(filepath), '\\.[^.]+$')
  tbl = read_tsv(filepath)
  tbl$Genome_id = filename
  return(tbl)
}

# Unique function with drop of NA and dash
unique_drop = function(values){
  values = unique(values)
  if( all(is.na(values)) ){
    return(NA)
  }else if(all(values %in% c('-', NA))){
    return('-')
  }else{
    values = values[!is.na(values)]
    values = values[values != '-']
    return(values)
  }
}

option_list = list(
  make_option(c("-f", "--fr_dict"), type = "character", 
              help = "File with functional roles dictionary."),
  make_option(c("-p", "--p_dict"), type = "character", 
              help = "File with phenotypes dictionary."),
  make_option(c("-a", "--ann_dir"), type = "character", 
              help = "Directory with annotations."),
  make_option(c("-o", "--out_dir"), type = "character",
              default = "./annotation_simplified",
              help = "Output_directory. Default: './annotation_simplified'.")
)

options = parse_args(OptionParser(option_list=option_list))

# options = parse_args(OptionParser(option_list=option_list),
#                      args = c("-f", "./dictionaries/functional_roles.txt",
#                               "-p", "./dictionaries/phenotypes.txt",
#                               "-a", "annotation"))

# Create output directory
dir.create(options$out_dir, showWarnings = FALSE)

# Read and concatenate annotation files
ann_files = Sys.glob(file.path(options$ann_dir, "*"))
ann_table = ann_files %>% map_dfr(read_tsv_name)

if(ann_table %>% nrow() == 0 ){
  print("No annotations were loaded. Please check path to ypur annotation directory.")
  quit(save = "no")
}

ann_table = ann_table %>% select(ID, Winner, Genome_id) %>% distinct()

# Read functional roles dictionary
fr_dict = read_tsv(options$fr_dict) %>% select(Role, Name, Phenotype)

# Read phenotypes dictionary
p_dict = read_tsv(options$p_dict) %>% select(Category, Description, Phenotype)

# Combine tables
combined_table = left_join(ann_table, fr_dict, by = c("Winner" = "Role"))

# Return Table with NA phenotypes
combined_table %>% filter(is.na(Phenotype)) %>% 
  rename(Gene_id=ID, Role=Winner) %>%
  select(Role) %>% 
  distinct() %>%
  write_tsv(file.path(options$out_dir, "no_phenotype.tsv"))

combined_table = combined_table %>% 
  mutate(Phenotype = str_split(Phenotype, '; ')) %>% 
  unnest(Phenotype)

combined_table = left_join(combined_table, p_dict, by = "Phenotype")

# Return Table with NA category
combined_table %>% filter(is.na(Category)) %>% 
  select(Phenotype) %>% 
  filter(!is.na(Phenotype), Phenotype != '-') %>%
  distinct() %>%
  write_tsv(file.path(options$out_dir, "no_category.tsv"))

# Collapse rows using "; " separator
combined_table = combined_table %>% 
  group_by(ID, Winner, Genome_id) %>%
  summarize(
    Name = paste0(unique(Name), collapse="; "),
    Phenotype = paste0(unique_drop(Phenotype), collapse="; "),
    Category = paste0(unique_drop(Category), collapse="; "),
    Description = paste0(unique_drop(Description), collapse="; ")
  )

# Set column names
combined_table = combined_table %>% rename(Gene_id=ID, Role=Winner) %>%
  select(Genome_id, Gene_id, Name, Role, Category, Description, Phenotype)

# Write output
combined_table %>% write_tsv(file.path(options$out_dir, "annotations_simple.tsv"))

