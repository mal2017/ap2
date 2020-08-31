cat PRJNA392905.yaml | grep "ftp" | cut -f 2 -d "[" | cut -f 1 -d ']' | grep ftp | while read FILE; do wget -P /scratch/mal456/fastq/immgen/ $FILE; done
