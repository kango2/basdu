#!/bin/bash

#resources: 8 threads
#input: bamfile, sampleid, SMRTCell ID, Well ID, rundate, sequencing centre, output directory


usage() { echo -e "\nUsage: $0 -h help -t <threads> -i <inputbam> -s <sampleid> -c <smrtcellid> -w <wellid> -o <outputdir> -d <rundate> -f <sequencingcentre>\n\n" 1>&2; exit 1; }
no_args="true"
while getopts ":hi:s:c:w:o:t:d:f:" option; do
    case "${option}" in
				h) usage;;
        i) inputbam=${OPTARG};;
        s) sampleid=${OPTARG};;
        c) smrtcellid=${OPTARG};;
        w) wellid=${OPTARG};;
        o) outputdir=${OPTARG};;
				t) threads=${OPTARG};;
				f) seqcentre=${OPTARG};;
				d) rundate=${OPTARG};;
				:) printf "missing argument for -%s\n" "$OPTARG" >&2; usage;;
			 \?) printf "illegal option: -%s\n" "$OPTARG" >&2; usage;;
        *) usage;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && usage
[[ -f "${inputbam}" ]] || usage
[[ -z "${sampleid}" ]] && usage
[[ -z "${smrtcellid}" ]] && usage
[[ -z "${wellid}" ]] && usage
[[ -z "${outputdir}" ]] && usage
[[ -z "${seqcentre}" ]] && usage
if [[ ! -d "${outputdir}" ]]
then
	mkdir -p "${outputdir}"
fi
if [[ -z "${threads}" ]]; then
	threads=1
fi
re='^[0-9]{1,2}$'
if ! [[ "${threads}" =~ $re ]] ; then
   echo "error: Not a number :$threads:" >&2;
	 usage
fi
re='^20[0-9]{6}$'
if ! [[ ${rundate} =~ $re ]] ; then
   echo "error: Not a valid date. Supply in YYYYMMDD format" >&2;
	 usage
fi
if [[ ${threads} -gt 48 ]]
then
	echo "${threads} set higher than 48. Reducing to 48 only."
	threads=48
fi

set -e
set -o pipefail

module load samtools cutadapt

chkpoint="${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.bam2fastx.done"
running="${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.bam2fastx.running"

if [[ -e "${running}" ]]
then
	exit 0
fi

if [[ ! -e "${chkpoint}"  || ! -s "${chkpoint}" || $(tail -n1 ${chkpoint} | cut -f3 -d',') != " EXIT_STATUS:0" ]]
then
	touch $running
	echo JOBID:"$JOBID", STARTTIME:$(date +%s) >> ${chkpoint}
	stime=$(date +%s)

samtools fastq --threads ${threads} ${inputbam} |\
cutadapt --cores ${threads} --anywhere file:${PACBIOADAPTERS} \
--error-rate 0.1 --overlap 25 --match-read-wildcards --revcomp --discard-trimmed \
--json ${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.cutadapt.json - 2>${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.cutadapt.txt |\
tee >(md5sum | sed "s/-/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.fq.data/" >>${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.md5) \
>(pigz -p ${threads} >${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.fq.gz) |\
sed -n '1~4s/^@/>/p;2~4p' | tee >(md5sum | sed "s/-/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.fa.data/" >>${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.md5) |\
pigz -p ${threads} >${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.fa.gz

md5sum ${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.fq.gz >>${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.md5
md5sum ${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.fa.gz >>${outputdir}/${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.md5

  exitstatus=$?
	etime=$(date +%s)
  echo JOBID:${PBS_JOBID}, TASKNAME:${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.bam2fastx, EXIT_STATUS:${exitstatus}, STARTTIME:${stime}, ENDTIME:${etime}, ELAPSED:$((etime - stime)) >>${chkpoint}

	rm -f $running
  if [ "${exitstatus}" -ne 0 ]
  then
    echo ERROR: ${sampleid}-${smrtcellid}-${wellid}-${seqcentre}-${rundate}.bam2fastx failed with ${exitstatus}
    exit ${exitstatus}
  fi
fi
