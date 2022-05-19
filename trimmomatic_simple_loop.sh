for infile in *1_001.fastq.gz
 do
   base=$(basename ${infile} 1_001.fastq.gz)
   java -jar /programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${infile} ${base}2_001.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 \
		-threads 10
 done

