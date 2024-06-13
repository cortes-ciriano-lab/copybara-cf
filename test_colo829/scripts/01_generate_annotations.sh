#!/bin/bash

##################################################
# Script to generate reference bin annotations 
# using different bin sizes
# Testing COPYBARA on ONT Colo829 matched normal data downsampled to 10x

### 2024/06/12
### csauer
##################################################


workdir="/nfs/research/icortes/csauer/projects/cfdna_ont_copy_number/COPYBARA/test_colo829"
scripts="${workdir}/../COPYBARA"

#in
ref="/nfs/research/icortes/DATA/hg38/Homo_sapiens_assembly38.fasta"
blacklist="/nfs/research/icortes/csauer/projects/QDNAseq.hg38_build/scripts/resources/ENCFF356LFX_hg38_blacklist.bed"

#out
annot_dir="${workdir}/ref_annotations"
logs="${annot_dir}/logs"

mkdir -p $annot_dir $logs

# loop through different window sizes...
for bin in 1 10 50 100 500; do
        echo "window size ${bin}kbp"

        annot="${annot_dir}/${bin}kbp_bin_ref_all.bed"

	# submit python scripts
        sbatch --mem=10G -t 12:00:00 -c 24 -n 1 -J annot_${bin}kbp -o ${logs}/annot_${bin}kbp.%J.o.log -e ${logs}/annot_${bin}kbp.%J.e.log << EOF
#!/bin/bash
date
python ${scripts}/bin_ref_generator_v3.py -w ${bin} -f ${ref} -b ${blacklist} -o ${annot_dir}
date
EOF
done
