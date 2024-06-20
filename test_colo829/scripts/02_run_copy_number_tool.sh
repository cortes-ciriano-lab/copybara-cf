#!/bin/bash

##################################################
# Script to run/test copy number tool on tumour ONT bams
# prelim analysis for Savana

### 2024/05/21
### csauer
##################################################

workdir="/nfs/research/icortes/csauer/projects/cfdna_ont_copy_number/COPYBARA/test_colo829"
scripts="${workdir}/../COPYBARA"

#in
nbam="/nfs/research/icortes/csauer/projects/cfdna_ont_copy_number/colo829_ont_test_data/input_data/downsampled_10x_bams_orig/colo829_normal_ont_10x_sorted.bam"
tbam="/nfs/research/icortes/csauer/projects/cfdna_ont_copy_number/colo829_ont_test_data/input_data/downsampled_10x_bams_orig/colo829_tumour_ont_10x_sorted.bam"
sn="colo829_10x_ont"
nmode="mnorm"

annot_dir="${workdir}/ref_annotations"

#out
outdir="${workdir}/copy_number_out_v3"
logs="${outdir}/logs"

mkdir -p $outdir $logs

IFS=$'\n'

# loop through different window sizes...
for bin in 500 100 50 10 1; do
	echo "window size ${bin}kbp"

	annot="${annot_dir}/${bin}kbp_bin_ref_all.bed"

	out="${outdir}/${bin}kbp"
	sample="${sn}_${bin}kbp"
	echo ${sample}

	mkdir -p ${out}
		
	# submit python scripts
	sbatch --mem=30G -t 96:00:00 -c 24 -J CN_${sample} -o ${logs}/CN_${sample}.%J.o.log -e ${logs}/CN_${sample}.%J.e.log << EOF
#!/bin/bash
date
python ${scripts}/bin_read_counter_v3.py -b ${tbam} -nb ${nbam} -s ${sample} -a ${annot} -o ${out} -t 24

python ${scripts}/smooth_CNdata_v2.py -rc ${out}/${sample}_read_counts_${nmode}_log2r.tsv -o ${out}

python ${scripts}/segment_CNdata_v2.py -i ${out}/${sample}_read_counts_${nmode}_log2r_smoothened_sl10_t0.025.tsv -o ${out} -t 24 -s 10000

date
EOF
done



