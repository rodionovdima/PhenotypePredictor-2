library(tidyverse)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-d", "--data"), action="store", default='./data', type='character', help="Path to the data directory with functional roles and phenotype dictionaries. Default: ./data")
parser <- add_option(parser, c("-a", "--annotation"), action="store", default='./annotation', type='character', help="Path to the directory with Danatello annotations. Default: ./annotation")
parser <- add_option(parser, c("-o", "--out"), action="store", default='./output', type='character', help="Output folder. Default: ./output")
args = parse_args(parser)

# Test:
# args = list(data = './data', annotations = './test/annotation', out = './test/out')

# Load Functional roles
print('Load Functional roles')
fr <- read_tsv(file.path(args$data, 'fr.txt'), show_col_types = F)

# Make name dictionary and select first name to unify names
fr_names <- fr %>% select(Name, Role) %>% 
  group_by(Role) %>% 
  mutate(N = row_number()) %>% 
  filter(N == 1) %>% 
  select(-N)

# Rename functional roles and remove duplicates and subsystems
fr <- fr %>% select(-Name, -Subsystem) %>% 
  left_join(fr_names, by = 'Role') %>% 
  distinct()

# Load phenotypes data
print('Load phenotypes data')
phenotype <- read_tsv(file.path(args$data, 'phenotype.txt'), show_col_types = F)

# Load annotations
read_annotation <- function(ann_file){
  print(paste('Reading file:', ann_file))
  ann_tbl <- read_tsv(ann_file, show_col_types = F)
  ann_tbl$Genome_name <- basename(ann_file)
  return(ann_tbl)
}

print('Load annotations')
ann <- Sys.glob(file.path(args$annotation, '*')) %>% map_dfr(read_annotation)

# Check that annotations exist
if( nrow(ann) == 0 ){
  stop("No annotations found.")
}

# Add genome names
if( file.exists(file.path(args$data, 'genome_names.txt')) ){
  print("Substitute file names with genome names")
  genome_table <- read_tsv(file.path(args$data, 'genome_names.txt'), col_types = cols(.default = col_character()))
  if( nrow(genome_table) > 0){
    absent_genomes <- setdiff(
      unique(str_remove(ann$Genome_name, '\\.fa.*$')),
      genome_table$Genome_ID
    )
    
    if( length(absent_genomes) > 0 ){
      stop(cat("Some genomes are not present in the genome_names.txt file:", absent_genomes, collapse = '\n'))
    }
    
    ann <- ann %>% mutate(Genome_name = str_remove(Genome_name, '\\.fa.*$')) %>% 
      rename(Genome_ID = Genome_name)
    ann <- ann %>% left_join(genome_table, by = "Genome_ID") %>% 
      select(-Genome_ID)
  }else{
    stop(paste("File", file.path(args$data, 'genome_table.txt'), "exists but empty"))
  }
}else{
  warning("No genome names were given. Leaving filenames.")
}

# Make output directory
dir.create(args$out, showWarnings = F)

# Write output file
print('Writing output file')
ann %>% 
  rename(Role = Winner) %>% 
  select(Genome_name, ID, Role) %>% 
  mutate(Role = str_remove(Role, ' #.+')) %>%
  mutate(Role = str_replace(Role, ' @ ', ' / ')) %>%
  mutate(Role = str_split(Role, ' / ')) %>% 
  unnest(Role) %>%
  distinct() %>% 
  left_join(fr, by = 'Role') %>% 
  left_join(phenotype, by = 'Role') %>% 
  relocate(Name, .after = 'ID') %>% 
  filter(!is.na(Module1)) %>% 
  write_tsv(file.path(args$out, 'out.txt'))
