# dg_merge_niftis - code to merge different images based on clusters saved after 2nd-level analysis in SPM.
# I modified each cluster/nifti with the modified # version of mw_tidy, giving voxels of every cluster a 
# different value (1 for c1,2 for c2,3 for c3,and so on). Use merged nifti with labels saved in a txt file
# for further FC analysis in CONN, following https://www.youtube.com/watch?v=Wjp_GOFD3oA&ab_channel=AndrewJahn

# pip install nibabel # WHEN NEEDED!

import nibabel as nib # https://nipy.org/nibabel/
import numpy as np

# list of NIfTI files
input_files = ['cluster1_both_AF_FWE_c0_h0_t0_s0_b1.nii','cluster2_both_AF_FWE_c0_h0_t0_s0_b1.nii','cluster3_both_AF_FWE_c0_h0_t0_s0_b1.nii', 
                'cluster4_both_AF_FWE_c0_h0_t0_s0_b1.nii','cluster5_both_AF_FWE_c0_h0_t0_s0_b1.nii','cluster6_both_AF_FWE_c0_h0_t0_s0_b1.nii',
                'cluster7_both_AF_FWE_c0_h0_t0_s0_b1.nii','cluster8_both_AF_FWE_c0_h0_t0_s0_b1.nii','cluster9_both_AF_FWE_c0_h0_t0_s0_b1.nii',
                'cluster10_both_AF_FWE_c0_h0_t0_s0_b1.nii']

# load the first file to get the shape and affine
first_img = nib.load(input_files[0])
merged_data = np.zeros(first_img.shape)

# loop through each file and add its data to the merged array
for i, file in enumerate(input_files):
    img = nib.load(file)
    data = img.get_fdata()

# assign the value (i+1) to the corresponding voxels
merged_data[data > 0] = i + 1

# create a new nifti with the merged data/clusters
merged_img = nib.Nifti1Image(merged_data, affine=first_img.affine)

# save it
nib.save(merged_img, 'merged_clusters.nii')
print("Merged NIfTI file saved as 'merged_clusters.nii'")
