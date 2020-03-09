library(tidyverse)
library(yaml)
library(readxl)

manifest <- read_tsv("config/PRJNA392905.txt") %>%
  dplyr::select(sample_title,fastq_ftp) %>%
  separate(fastq_ftp,into=c("r1","r2"),sep=";") %>%
  mutate(sample_title = str_replace_all(sample_title,"-","_")) %>%
  mutate(SampleName = sample_title) %>%
  separate(sample_title,into = c("biosample","replicate"), sep="\\#")

read_excel("config/NIHMS1538358-supplement-8.xlsx") %>%
  dplyr::select(SampleName, CellType, ImmGenLab, Lineage, CellFamily, Organ) %>%
  left_join(manifest,.) %>%
  drop_na() %>%
  group_by(SampleName) %>%
  summarize(fastq= pmap(list(r1,r2,CellType,Lineage,CellFamily,Organ), .f= function(x,y,z,a,b,c) list(fastq=list(r1=x,r2=y),condition=list(z,a,b,c)))) %>%
  mutate(SampleName=str_replace(SampleName,"\\#","\\.")) %>%
  deframe() %>%
  yaml::write_yaml("config/PRJNA392905.yaml")
