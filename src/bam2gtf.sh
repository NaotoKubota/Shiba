#!/bin/bash -e
set -o pipefail
# Prevent commands misbehaving due to locale differences
export LC_ALL=C

function usage {
    cat <<EOS

Usage: $(basename "$0") -i experiment.txt -r reference_annotation.gtf -o assembled_annotation.gtf -p [VALUE]

    -h  Display help
    -i  Experiment table
    -r  Reference GTF file
    -o  Assembled GTF file
    -p	Number of processors to use (default: 1)

EOS
    exit 2
}


PROCESS=1


function lack_of_necessary_param() {
    usage
    exit 1
}


IS_THERE_NECESSARY_OPT_i=false
IS_THERE_NECESSARY_OPT_r=false
IS_THERE_NECESSARY_OPT_o=false


while getopts "i:r:o:p:h" optKey; do
    case "$optKey" in
    i)
        IS_THERE_NECESSARY_OPT_i=true
        EXPERIMENT=${OPTARG}
        ;;
    r)
        IS_THERE_NECESSARY_OPT_r=true
        REFERENCE=${OPTARG}
        ;;
    o)
        IS_THERE_NECESSARY_OPT_o=true
        ASSEMBLED=${OPTARG}
        ;;
	p)
        PROCESS=${OPTARG}
        ;;
    h|* )
        usage
        ;;
    esac
done


if [ "${IS_THERE_NECESSARY_OPT_i}" == false ] || [ "${IS_THERE_NECESSARY_OPT_r}" == false ] || [ "${IS_THERE_NECESSARY_OPT_o}" == false ]; then
    lack_of_necessary_param
fi;


echo -e "Start transcirpt assembly!"


OUTPUTDIR=$(echo ${ASSEMBLED} | xargs -n1 dirname)
mkdir -p ${OUTPUTDIR}
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

    # StringTie2
    stringtie \
    -p ${PROCESS} \
    -G ${REFERENCE} \
    -o ${DIR}/${SAMPLE}.assembled.gtf \
    ${BAM} && \
    echo ${DIR}/${SAMPLE}.assembled.gtf \
    >> ${OUTPUTDIR}/gtf_list.txt

done

echo -e "Merging each gtf..."

# Merge GTF
GTFLIST=$(cat ${OUTPUTDIR}/gtf_list.txt) && \
stringtie --merge \
-p ${PROCESS} \
-G ${REFERENCE} \
-o ${ASSEMBLED} \
${GTFLIST}

rm -rf ${OUTPUTDIR}/gtf_list.txt

echo -e "Done!"
