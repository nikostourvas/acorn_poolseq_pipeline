#!/bin/bash

JAVA_DIR=/usr/share/java
BAM=${1}
DEDUP_BAM=${2}

# remove duplicates with picard MarkDuplicates
java -jar -Xms256m -Xmx27g ${JAVA_DIR}/picard.jar MarkDuplicates \
	I=${BAM} O=${DEDUP_BAM} M=${DEDUP_BAM}.picardstat \
	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
	REMOVE_DUPLICATES=true

# gather statistics with samtools flagstat    
samtools flagstat ${DEDUP_BAM} > ${DEDUP_BAM}.flagstat