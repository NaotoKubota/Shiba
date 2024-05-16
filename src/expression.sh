#!/bin/bash -e
set -o pipefail
# Prevent commands misbehaving due to locale differences
export LC_ALL=C

function usage {
    cat <<EOS

Usage: $(basename "$0") -i experiment.txt -g reference_annotation.gtf -o output_dir -r [VALUE] -a [VALUE] -p [VALUE]

    -h  Display help
    -i  Experiment table
    -g  Reference GTF file
    -o 	Output directory
    -r  Reference group for differential expression analysis (default: NA)
    -a  Alternative group for differential expression analysis (default: NA)
    -p  Number of processors to use (default: 1)

EOS
    exit 2
}


SRC_PATH=$(dirname "$0")
PROCESS=1
REFGROUP="NA"
EXPGROUP="NA"


function lack_of_necessary_param() {
    usage
    exit 1
}


IS_THERE_NECESSARY_OPT_i=false
IS_THERE_NECESSARY_OPT_g=false
IS_THERE_NECESSARY_OPT_o=false


while getopts "i:g:o:r:e:p:h" optKey; do
    case "$optKey" in
    i)
        IS_THERE_NECESSARY_OPT_i=true
        EXPERIMENT=${OPTARG}
        ;;
    g)
        IS_THERE_NECESSARY_OPT_g=true
        REFERENCE=${OPTARG}
        ;;
    o)
        IS_THERE_NECESSARY_OPT_o=true
        OUTPUTDIR=${OPTARG}
        ;;
	r)
        REFGROUP=${OPTARG}
        ;;
    e)
        EXPGROUP=${OPTARG}
        ;;
    p)
        PROCESS=${OPTARG}
        ;;
    h|* )
        usage
        ;;
    esac
done


if [ "${IS_THERE_NECESSARY_OPT_i}" == false ] || [ "${IS_THERE_NECESSARY_OPT_g}" == false ] || [ "${IS_THERE_NECESSARY_OPT_o}" == false ]; then
    lack_of_necessary_param
fi;


mkdir -p ${OUTPUTDIR}/logs
rm -rf ${OUTPUTDIR}/logs/featureCounts.log
rm -rf ${OUTPUTDIR}/countfiles.txt
cat ${EXPERIMENT} | while read line || [ -n "${line}" ]
do

	SAMPLE=$(echo ${line} | cut -d " " -f 1)
	BAM=$(echo ${line} | cut -d " " -f 2)
	DIR=$(echo ${line} | cut -d " " -f 2 | xargs -n1 dirname)

    if [ "${SAMPLE}" == "sample" ]; then

        continue

    fi

    echo -e "${SAMPLE}...."

    # Check whether bam index file exists
    BAM_INDEX=${BAM}.bai
    if [ ! -f "${BAM_INDEX}" ]; then

        echo -e "Make bam index..." && \
        samtools index ${BAM}

    else

        :

    fi

	# Check if BAM is paired/single-end
	FLAG_PAIRED=$({ samtools view -H ${BAM}; samtools view ${BAM} | head -n 1; } | samtools view -c -f 1)

	# featureCounts
	if [[ ${FLAG_PAIRED} == 1 ]]; then

		# Paired-end
		featureCounts \
		-a ${REFERENCE} \
		-o ${OUTPUTDIR}/${SAMPLE}_counts.txt \
		-T ${PROCESS} \
		-t exon \
		-g gene_id \
		-p \
		-B \
		${BAM} 2>> ${OUTPUTDIR}/logs/featureCounts.log

	else

		# Single-end
		featureCounts \
		-a ${REFERENCE} \
		-o ${OUTPUTDIR}/${SAMPLE}_counts.txt \
		-T ${PROCESS} \
		-t exon \
		-g gene_id \
		${BAM} 2>> ${OUTPUTDIR}/logs/featureCounts.log

    fi

	cat ${OUTPUTDIR}/${SAMPLE}_counts.txt | \
	sed -e 1d | \
	awk -F"\t" -v OFS="\t" -v SAMPLE=${SAMPLE} 'NR == 1{print $1,$6,SAMPLE}NR != 1{print $1,$6,$7}' \
	> ${OUTPUTDIR}/${SAMPLE}_counts_simple.txt

	echo -e "${OUTPUTDIR}/${SAMPLE}_counts_simple.txt" >> ${OUTPUTDIR}/countfiles.txt

done


# Calculate TPM
python ${SRC_PATH}/tpm.py \
${OUTPUTDIR}/countfiles.txt \
${OUTPUTDIR}/


if [ "${REFGROUP}" == "NA" ] || [ "${EXPGROUP}" == "NA" ]; then

    :

else

    echo -e "Differential expression analysis by DESeq2..."
    # Differential expression analysis by DESeq2
    Rscript ${SRC_PATH}/deseq2.R \
    ${EXPERIMENT} \
    ${OUTPUTDIR}/counts.txt \
    ${REFGROUP} \
    ${EXPGROUP} \
    ${OUTPUTDIR}/DEG.txt 2> ${OUTPUTDIR}/logs/DESeq2.log

fi

rm -rf ${OUTPUTDIR}/countfiles.txt ${OUTPUTDIR}/*_counts_simple.txt ${OUTPUTDIR}/*_counts.txt*
