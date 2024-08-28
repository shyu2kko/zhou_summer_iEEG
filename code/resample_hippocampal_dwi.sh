#!/bin/bash

export ANTSPATH=/data/mica1/01_programs/ants-2.3.4/bin/
export ZBTOOLDIR=/host/verges/tank/data/yigu/z-brains/

# settint up
export outdir=/host/verges/tank/data/yigu/missing_hippocampi/
IFS=$'\n' read -d '' -r -a list_to_proc < /data/mica1/03_projects/yigu/summer_iEEG/data/everyone.txt

sesid="01"

for subjid in ${list_to_proc[@]}; do

	${ZBTOOLDIR}/zbrains \
		--run proc \
		--sub sub-${subjid} \
		--ses ses-${sesid} \
		--dataset /data/mica3/BIDS_MICs/ \
		--zbrains $outdir/zbrains \
		--micapipe /micapipe_v0.2.0/ \
		--hippunfold /hippunfold_v1.3.0/ \
		--struct hippocampus \
		--feat "ADC" "FA" "qT1"

done
