#!/bin/bash

main_files=(
"merged_rep_MDEC_est_5_JH_0345_E0419.dedup.overlapclipped.bam"
"merged_rep_MPHA_est_20_JH_0356_E0616.dedup.overlapclipped.bam"
"merged_rep_MPHA_est_20_JH_0360_E0623.dedup.overlapclipped.bam"
"merged_rep_MSEN_est_20_JH_0363_E0468.dedup.overlapclipped.bam"
"merged_rep_MPHA_est_10_JH_0284_E0413.dedup.overlapclipped.bam"
"merged_rep_MPHA_kal_20_JH_0215_E0674.dedup.overlapclipped.bam"
"merged_rep_MPHA_kal_20_JH_0214_E0396.dedup.overlapclipped.bam"
"merged_rep_MPHA_kal_20_JH_0211_E0246.dedup.overlapclipped.bam"
"merged_rep_MPHA_kal_20_JH_0206_E0296.dedup.overlapclipped.bam"
"merged_rep_MSEN_kal_10_JH_0085_E0667.dedup.overlapclipped.bam"
"merged_rep_MSEN_sna_20_JH_0052_E0357.dedup.overlapclipped.bam"
"merged_rep_MDEC_sna_20_JH_0042_E0280.dedup.overlapclipped.bam"
"merged_rep_MSEN_sna_20_JH_0041_E0584.dedup.overlapclipped.bam"
"merged_rep_MPHA_sna_5_JH_0594_E0718.dedup.overlapclipped.bam"
"merged_rep_MDEC_kal_10_JH_0130_E0716.dedup.overlapclipped.bam"
"merged_rep_MPHA_sna_20_JH_0715_E0717.dedup.overlapclipped.bam"
"merged_rep_MPHA_est_20_JH_0408_E0619.dedup.overlapclipped.bam"
"merged_rep_MDEC_kal_5_JH_0560_E0603.dedup.overlapclipped.bam"
"merged_rep_MCAR_kal_10_JH_0143_E0274.dedup.overlapclipped.bam"
"merged_rep_MSEN_kal_10_JH_0550_E0349.dedup.overlapclipped.bam"
"merged_rep_MMIR_sna_10_DV_0497_E0734.dedup.overlapclipped.bam"
"merged_clone_MDEC_kal_20_JH_0167_E0706.dedup.overlapclipped.bam"
"merged_clone_MDEC_kal_10_JH_0252_E0557.dedup.overlapclipped.bam"
"merged_clone_MDEC_est_10_JH_0390_E0342.dedup.overlapclipped.bam"
"merged_clone_MDEC_sna_5_JH_0604_E0500.dedup.overlapclipped.bam"
"merged_clone_MDEC_dir_5_JH_0660_E0446.dedup.overlapclipped.bam"
"merged_clone_MMIR_sna_10_DV_0337_E0757.dedup.overlapclipped.bam"
"merged_clone_MMIR_sna_10_DV_0338_E0737.dedup.overlapclipped.bam"
"merged_clone_MMIR_sna_10_DV_0424_E0762.dedup.overlapclipped.bam"
"merged_clone_MMIR_sna_20_DV_0275_E0768.dedup.overlapclipped.bam"
"merged_clone_MMIR_sna_20_DV_0451_E0740.dedup.overlapclipped.bam"
"merged_clone_MPHA_kal_20_JH_0193_E0338.dedup.overlapclipped.bam"
"merged_clone_MPHA_kal_20_JH_0222_E0494.dedup.overlapclipped.bam"
"merged_clone_MPHA_kal_10_JH_0262_E0243.dedup.overlapclipped.bam"
"merged_clone_MPHA_est_10_JH_0384_E0529.dedup.overlapclipped.bam"
"merged_clone_MPHA_est_20_JH_0467_E0709.dedup.overlapclipped.bam"
"merged_clone_MPHA_est_20_JH_0477_E0527.dedup.overlapclipped.bam"
"merged_clone_MPHA_est_20_JH_0484_E0653.dedup.overlapclipped.bam"
"merged_clone_MPHA_sea_20_JH_0516_E0650.dedup.overlapclipped.bam"
"merged_clone_MPHA_sna_10_JH_0589_E0590.dedup.overlapclipped.bam"
"merged_clone_MPHA_sna_5_JH_0596_E0415.dedup.overlapclipped.bam"
"merged_clone_MPHA_wat_20_JH_0619_E0104.dedup.overlapclipped.bam"
"merged_clone_MPHA_dir_20_JH_0652_E0130.dedup.overlapclipped.bam"
"merged_clone_MPHA_dir_20_JH_0658_E0126.dedup.overlapclipped.bam"
"merged_clone_MPHA_hul_20_JH_0678_E0559.dedup.overlapclipped.bam"
"merged_clone_MPHA_sna_20_JH_0695_E0277.dedup.overlapclipped.bam"
"merged_clone_MPHA_kal_20_JH_0723_E0678.dedup.overlapclipped.bam"
"merged_clone_MSEN_sna_10_JH_0012_E0470.dedup.overlapclipped.bam"
"merged_clone_MSEN_sna_10_JH_0031_E0621.dedup.overlapclipped.bam"
"merged_clone_MSEN_sna_10_JH_0011_E0361.dedup.overlapclipped.bam"
"merged_clone_MSEN_kal_10_JH_0124_E0066.dedup.overlapclipped.bam"
"merged_clone_MSEN_kal_20_JH_0213_E0017.dedup.overlapclipped.bam"
"merged_clone_MSEN_est_20_JH_0307_E0090.dedup.overlapclipped.bam"
"merged_clone_MSEN_est_20_JH_0462_E0399.dedup.overlapclipped.bam"
"merged_clone_MSEN_sna_20_JH_0709_E0100.dedup.overlapclipped.bam"

)

intermediate_files=(
"merged_MPHA_est_5_JH_0345_E0419.dedup.overlapclipped.bam merged_MDEC_est_5_JH_0420_E0113.dedup.overlapclipped.bam"
"merged_MPHA_est_20_JH_0404_E0088.dedup.overlapclipped.bam merged_MPHA_est_20_JH_0356_E0616.dedup.overlapclipped.bam"
"merged_MPHA_est_20_JH_0405_E0325.dedup.overlapclipped.bam merged_MPHA_est_20_JH_0360_E0623.dedup.overlapclipped.bam"
"merged_MSEN_est_20_JH_0409_E0573.dedup.overlapclipped.bam merged_MSEN_est_20_JH_0363_E0468.dedup.overlapclipped.bam"
"merged_MPHA_est_10_JH_0287_E0109.dedup.overlapclipped.bam merged_MPHA_est_10_JH_0284_E0413.dedup.overlapclipped.bam"
"merged_MPHA_kal_20_JH_0540_E0403.dedup.overlapclipped.bam merged_MPHA_kal_20_JH_0215_E0674.dedup.overlapclipped.bam"
"merged_MPHA_kal_20_JH_0539_E0068.dedup.overlapclipped.bam merged_MPHA_kal_20_JH_0214_E0396.dedup.overlapclipped.bam"
"merged_MPHA_kal_20_JH_0211_E0246.dedup.overlapclipped.bam merged_MPHA_kal_20_JH_0544_E0257.dedup.overlapclipped.bam"
"merged_MPHA_kal_20_JH_0545_E0197.dedup.overlapclipped.bam merged_MPHA_kal_20_JH_0206_E0296.dedup.overlapclipped.bam"
"merged_MSEN_kal_10_JH_0555_E0295.dedup.overlapclipped.bam merged_MSEN_kal_10_JH_0085_E0667.dedup.overlapclipped.bam"
"merged_MSEN_sna_20_JH_0052_E0357.dedup.overlapclipped.bam merged_MSEN_sna_20_JH_0700_E0535.dedup.overlapclipped.bam"
"merged_MDEC_sna_20_JH_0042_E0280.dedup.overlapclipped.bam merged_MDEC_sna_20_JH_0582_E0294.dedup.overlapclipped.bam"
"merged_MSEN_sna_20_JH_0714_E0475.dedup.overlapclipped.bam merged_MSEN_sna_20_JH_0041_E0584.dedup.overlapclipped.bam"
"merged_MDEC_sna_5_JH_0594_E0585.dedup.overlapclipped.bam merged_MDEC_sna_5_JH_0594_E0718.dedup.overlapclipped.bam"
"merged_MDEC_kal_10_JH_0130_E0644.dedup.overlapclipped.bam merged_MDEC_kal_10_JH_0130_E0716.dedup.overlapclipped.bam"
"merged_MPHA_sna_20_JH_0715_E0658.dedup.overlapclipped.bam merged_MPHA_sna_20_JH_0715_E0717.dedup.overlapclipped.bam"
"merged_MPHA_est_20_JH_0408_E0715.dedup.overlapclipped.bam merged_MPHA_est_20_JH_0408_E0619.dedup.overlapclipped.bam"
"merged_MDEC_kal_5_JH_0560_E0719.dedup.overlapclipped.bam merged_MDEC_kal_5_JH_0560_E0603.dedup.overlapclipped.bam"
"merged_MCAR_kal_10_JH_0143_E0274.dedup.overlapclipped.bam merged_MCAR_kal_10_JH_0136_E0433.dedup.overlapclipped.bam"
"merged_MSEN_kal_10_JH_0550_E0349.dedup.overlapclipped.bam merged_MDEC_kal_10_JH_0089_E0485.dedup.overlapclipped.bam"
"merged_MMIR_sna_20_DV_0370_E0761.dedup.overlapclipped.bam merged_MMIR_sna_20_DV_0452_E0741.dedup.overlapclipped.bam"
"merged_MMIR_sna_10_DV_0497_E0734.dedup.overlapclipped.bam merged_MMIR_sna_10_DV_0498_E0735.dedup.overlapclipped.bam"
"merged_UNKN_kal_20_JH_0167_E0706.dedup.overlapclipped.bam merged_MDEC_kal_20_JH_0724_E0310.dedup.overlapclipped.bam"
"merged_UNKN_kal_10_JH_0252_E0557.dedup.overlapclipped.bam merged_MDEC_kal_10_JH_0547_E0677.dedup.overlapclipped.bam"
"merged_MDEC_est_10_JH_0390_E0342.dedup.overlapclipped.bam merged_MDEC_est_10_JH_0387_E0033.dedup.overlapclipped.bam"
"merged_MDEC_sna_5_JH_0604_E0500.dedup.overlapclipped.bam merged_MDEC_sna_5_JH_0593_E0297.dedup.overlapclipped.bam"
"merged_MDEC_dir_5_JH_0660_E0446.dedup.overlapclipped.bam merged_MDEC_dir_5_JH_0663_E0477.dedup.overlapclipped.bam"
"merged_MMIR_sna_10_DV_0337_E0757.dedup.overlapclipped.bam merged_MMIR_sna_10_DV_0325_E0754.dedup.overlapclipped.bam"
"merged_MMIR_sna_10_DV_0338_E0737.dedup.overlapclipped.bam merged_MMIR_sna_10_DV_0331_E0766.dedup.overlapclipped.bam"
"merged_MMIR_sna_10_DV_0424_E0762.dedup.overlapclipped.bam merged_MMIR_sna_10_DV_0416_E0758.dedup.overlapclipped.bam"
"merged_MMIR_sna_20_DV_0275_E0768.dedup.overlapclipped.bam merged_MMIR_sna_20_DV_0271_E0760.dedup.overlapclipped.bam"
"merged_MMIR_sna_20_DV_0451_E0740.dedup.overlapclipped.bam merged_rep_MMIR_sna_20_DV_0370_E0761.dedup.overlapclipped.bam"
"merged_UNKN_kal_20_JH_0193_E0338.dedup.overlapclipped.bam merged_UNKN_kal_20_JH_0192_E0070.dedup.overlapclipped.bam"
"merged_UNKN_kal_20_JH_0222_E0494.dedup.overlapclipped.bam merged_UNKN_kal_20_JH_0221_E0307.dedup.overlapclipped.bam"
"merged_UNKN_kal_10_JH_0262_E0243.dedup.overlapclipped.bam merged_MPHA_kal_10_JH_0079_E0631.dedup.overlapclipped.bam"
"merged_MPHA_est_10_JH_0384_E0529.dedup.overlapclipped.bam merged_MPHA_est_10_JH_0381_E0505.dedup.overlapclipped.bam"
"merged_UNKN_est_20_JH_0467_E0709.dedup.overlapclipped.bam merged_UNKN_est_20_JH_0466_E0458.dedup.overlapclipped.bam"
"merged_UNKN_est_20_JH_0477_E0527.dedup.overlapclipped.bam merged_UNKN_est_20_JH_0476_E0360.dedup.overlapclipped.bam merged_UNKN_est_20_JH_0475_E0074.dedup.overlapclipped.bam"
"merged_UNKN_est_20_JH_0484_E0653.dedup.overlapclipped.bam merged_UNKN_est_20_JH_0483_E0423.dedup.overlapclipped.bam"
"merged_MPHA_sea_20_JH_0516_E0650.dedup.overlapclipped.bam merged_MDEC_sea_20_JH_0517_E0496.dedup.overlapclipped.bam"
"merged_MPHA_sna_10_JH_0589_E0590.dedup.overlapclipped.bam merged_MPHA_sna_10_JH_0062_E0686.dedup.overlapclipped.bam"
"merged_MSEN_sna_5_JH_0596_E0415.dedup.overlapclipped.bam merged_MPHA_sna_5_JH_0033_E0425.dedup.overlapclipped.bam"
"merged_MPHA_wat_20_JH_0619_E0104.dedup.overlapclipped.bam merged_MPHA_wat_20_JH_0617_E0565.dedup.overlapclipped.bam"
"merged_MPHA_dir_20_JH_0652_E0130.dedup.overlapclipped.bam merged_MPHA_dir_20_JH_0651_E0613.dedup.overlapclipped.bam"
"merged_MPHA_dir_20_JH_0658_E0126.dedup.overlapclipped.bam merged_MPHA_dir_20_JH_0657_E0422.dedup.overlapclipped.bam"
"merged_MCAR_hul_20_JH_0678_E0559.dedup.overlapclipped.bam merged_MCAR_hul_20_JH_0677_E0697.dedup.overlapclipped.bam"
"merged_MPHA_sna_20_JH_0695_E0277.dedup.overlapclipped.bam merged_MPHA_sna_20_JH_0694_E0546.dedup.overlapclipped.bam merged_MPHA_sna_20_JH_0696_E0502.dedup.overlapclipped.bam"
"merged_MPHA_kal_20_JH_0723_E0678.dedup.overlapclipped.bam merged_MPHA_kal_20_JH_0542_E0160.dedup.overlapclipped.bam"
"merged_MPHA_sna_10_JH_0012_E0470.dedup.overlapclipped.bam merged_SINT_sna_20_JH_0013_E0713.dedup.overlapclipped.bam"
"merged_MSEN_sna_10_JH_0031_E0621.dedup.overlapclipped.bam merged_MPHA_sna_20_JH_0707_E0201.dedup.overlapclipped.bam"
"merged_MPHA_sna_10_JH_0011_E0361.dedup.overlapclipped.bam merged_MPHA_sna_10_JH_0071_E0089.dedup.overlapclipped.bam"
"merged_MSEN_kal_10_JH_0124_E0066.dedup.overlapclipped.bam merged_MSEN_kal_10_JH_0122_E0515.dedup.overlapclipped.bam"
"merged_UNKN_kal_20_JH_0213_E0017.dedup.overlapclipped.bam merged_UNKN_kal_20_JH_0212_E0533.dedup.overlapclipped.bam"
"merged_UNKN_est_20_JH_0307_E0090.dedup.overlapclipped.bam merged_UNKN_est_20_JH_0306_E0261.dedup.overlapclipped.bam"
"merged_MPHA_est_20_JH_0462_E0399.dedup.overlapclipped.bam merged_UNKN_est_20_JH_0464_E0572.dedup.overlapclipped.bam"
"merged_MSEN_sna_20_JH_0709_E0100.dedup.overlapclipped.bam merged_MDEC_sna_20_JH_0049_E0428.dedup.overlapclipped.bam"
)

for ((i=0; i<${#main_files[@]}; i++)); do
  main_file="${main_files[i]}"
  intermediates="${intermediate_files[i]}"

  if [ -f "$main_file" ]; then
    echo "$main_file exists."
    for intermediate in $intermediates; do
      if [ -f "$intermediate" ]; then
        rm "$intermediate"
        echo "Deleted $intermediate."
      else
        echo "$intermediate does not exist."
      fi
    done
  else
    echo "$main_file does not exist."
  fi
done

# Remove potential clones that are slightly different than their clone group
rm merged_MPHA_wat_20_JH_0618_E0233.dedup.overlapclipped.bam
rm merged_UNKN_est_20_JH_0463_E0322.dedup.overlapclipped.bam
rm merged_MDEC_dir_5_JH_0661_E0708.dedup.overlapclipped.bam
rm merged_MMIR_sna_10_DV_0311_E0764.dedup.overlapclipped.bam
rm merged_MMIR_sna_10_DV_0323_E0753.dedup.overlapclipped.bam

# Remove 3 SINT
rm merged_UNKN_kal_40_JH_0147_E0236.dedup.overlapclipped.bam
rm merged_UNKN_kal_40_JH_0105_E0657.dedup.overlapclipped.bam
rm merged_UNKN_sna_60_JH_0075_E0215.dedup.overlapclipped.bam


