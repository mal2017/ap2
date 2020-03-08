library(tidyverse)
library(yaml)


manifest <- read_tsv("config/PRJNA301969.txt") %>%
  dplyr::select(sample_title,fastq_ftp) %>%
  separate(fastq_ftp,into=c("r1","r2"),sep=";") %>%
  mutate(sample_title = str_replace_all(sample_title,"-","_"))

read_yaml("config/heme-meta.yaml") %>% enframe() %>%
  mutate(condition = map_chr(value,`[[`,1)) %>%
  mutate(donor = map_chr(value,`[[`,2)) %>%
  dplyr::select(name, condition) %>%
  left_join(manifest,., by=c("sample_title"="name")) %>%
  drop_na() %>%
  group_by(sample_title) %>%
  summarize(fastq= pmap(list(r1,r2,condition), .f= function(x,y,z) list(fastq=list(r1=x,r2=y),condition=list(z,"primary")))) %>%
  deframe() %>%
  yaml::write_yaml("config/PRJNA301969.yaml")
