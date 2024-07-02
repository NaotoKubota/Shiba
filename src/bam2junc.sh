#!/bin/bash -e
set -o pipefail
# Prevent commands misbehaving due to locale differences
export LC_ALL=C

function usage {
    cat <<EOS

Usage: $(basename "$0") -i experiment.tsv -r RI_EVENT.txt -o junctions.bed -p [VALUE] -a [VALUE] -m [VALUE] -M [VALUE] -s [VALUE]

    -h  Display help
    -i  Experiment table
    -r  Intron retention event
    -o  Junction read counts
    -p  Number of processors to use (default: 1)
    -a  Minimum anchor length (default: 8)
    -m  Minimum intron size (default: 70)
    -M  Maximum intron size (default: 500000)
    -s  Strand specificity of RNA library preparation; XS: unstranded, RF: first-strand, FR: second-strand (default: XS)

EOS
    exit 2
}

SRC_PATH=$(dirname "$0")
ANCHOR=8
MINIMUM=70
MAXIMUM=500000
STRAND=XS
PROCESS=1

function lack_of_necessary_param() {
    usage
    exit 1
}

IS_THERE_NECESSARY_OPT_i=false
IS_THERE_NECESSARY_OPT_r=false
IS_THERE_NECESSARY_OPT_o=false

while getopts "i:a:m:M:s:r:o:p:h" optKey; do
    case "$optKey" in
    i)
        IS_THERE_NECESSARY_OPT_i=true
        EXPERIMENT=${OPTARG}
        ;;
    r)
        IS_THERE_NECESSARY_OPT_r=true
        RI_EVENT=${OPTARG}
        ;;
    o)
        IS_THERE_NECESSARY_OPT_o=true
        OUTPUT=${OPTARG}
        ;;
    p)
        PROCESS=${OPTARG}
        ;;
    a)
        ANCHOR=${OPTARG}
        ;;
    m)
        MINIMUM=${OPTARG}
        ;;
    M)
        MAXIMUM=${OPTARG}
        ;;
    s)
        STRAND=${OPTARG}
        ;;
    h|* )
        usage
        ;;
    esac
done


if [ "${IS_THERE_NECESSARY_OPT_i}" == false ] || [ "${IS_THERE_NECESSARY_OPT_r}" == false ] || [ "${IS_THERE_NECESSARY_OPT_o}" == false ]; then
    lack_of_necessary_param
fi;

OUTPUTDIR=$(dirname "${OUTPUT}")
mkdir -p ${OUTPUTDIR}/logs

cat ${RI_EVENT} | \
cut -f 6,7 | \
sed -e 1d | \
awk -F"\t" -v OFS="\t" '{split($1,l,":"); split(l[2],m,"-"); print l[1]":"m[1]"-"m[1]+1,l[1],m[1],m[1]+1,$2; print l[1]":"m[2]-1"-"m[2],l[1],m[2]-1,m[2],$2}' | \
awk '!a[$0]++' \
> ${OUTPUTDIR}/RI.saf

rm -rf ${OUTPUTDIR}/logs/featureCounts.log ${OUTPUTDIR}/logs/regtools.log
rm -rf ${OUTPUTDIR}/juncfiles.txt
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

    # exon-exon juntions
    # echo -e "Counting exon-exon junction reads by regtools...."

    regtools junctions extract \
    -s ${STRAND} \
    -a ${ANCHOR} \
    -m ${MINIMUM} \
    -M ${MAXIMUM} \
    -o ${DIR}/${SAMPLE}_exon-exon.junc \
    ${BAM} 2>> ${OUTPUTDIR}/logs/regtools.log

    echo -e "${DIR}/${SAMPLE}_exon-exon.junc\texon-exon" >> ${OUTPUTDIR}/juncfiles.txt

    # exon-intron junctions
    # echo -e "Counting exon-intron junction reads by featureCounts...."

    # Check if BAM is paired/single-end
    FLAG_PAIRED=$({ samtools view -H ${BAM}; samtools view ${BAM} | head -n 1; } | samtools view -c -f 1)

    if [[ ${FLAG_PAIRED} == 1 ]]; then

        featureCounts \
        -a ${OUTPUTDIR}/RI.saf \
        -o ${DIR}/${SAMPLE}_exon-intron.junc \
        -F SAF \
        --fracOverlapFeature 1.0 \
        -T ${PROCESS} \
        -O \
        -p \
        ${BAM} 2>> ${OUTPUTDIR}/logs/featureCounts.log

    else

        featureCounts \
        -a ${OUTPUTDIR}/RI.saf \
        -o ${DIR}/${SAMPLE}_exon-intron.junc \
        -F SAF \
        --fracOverlapFeature 1.0 \
        -T ${PROCESS} \
        -O \
        ${BAM} 2>> ${OUTPUTDIR}/logs/featureCounts.log

    fi

    echo -e "${DIR}/${SAMPLE}_exon-intron.junc\texon-intron" >> ${OUTPUTDIR}/juncfiles.txt

done

# Merge junction read counts table
python ${SRC_PATH}/merge_junc.py \
${OUTPUTDIR}/juncfiles.txt \
${OUTPUT} && \
rm -rf ${OUTPUTDIR}/juncfiles.txt ${OUTPUTDIR}/RI.saf
