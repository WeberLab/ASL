#!/usr/bin/env python3

import nibabel as nib
import numpy as np
import nipype.interfaces.fsl as fsl
import sys
import os.path as path

def make_cbf(fm0,fpcasl,PLD,SliceDelay,LabelDuration):
	#--------------------------------------------------------------------------
	## Parameters needed for the quantification. Most are not present in the
	## DICOM header
	#-------------------------------------------------------------------------
	Lambda=0.9 #Blood-brain partition coefficient.
	T1b=1.65 #T1 of the blood at 3T.
	alpha=0.85 # Labeling efficiency for a pCASL sequence.

	# # ACQUISITION PARAMS MUST BE CONFIRMED FOR EACH STUDY
	# # MINT PARAMS:
	# PLD=1.60 # Post-label delay (seconds) confirmed in ExamCard MJ
	# SliceDelay=0.039 #Slice delay. Confirmed on scanner per Guillaume's instructs MJ
	# LabelDuration=1.65 #Label duration. (confirmed in ExamCard MJ)

	# Load files
	print('Loading image files')

	m0 = nib.parrec.load(fm0,scaling='fp',strict_sort=True) # 3D Ref
	pcasl = nib.parrec.load(fpcasl,scaling='fp',strict_sort=True) # 4D pcasl

	print('Computing Scaling Factors')

	if len(set(pcasl.header.image_defs['scale slope'])) > 1:
		print('Warning! Multiple scale slope values for asl. Using first')
	if len(set(pcasl.header.image_defs['rescale slope'])) > 1:
		print('Warning! Multiple rescale slope values for asl. Using first')
	if len(set(pcasl.header.image_defs['rescale intercept'])) > 1:
		print('Warning! Multiple rescale intercept values for asl. Using first')
	if len(set(m0.header.image_defs['scale slope'])) > 1:
		print('Warning! Multiple scale slope values for M0. Using first')
	if len(set(m0.header.image_defs['rescale slope'])) > 1:
		print('Warning! Multiple rescale slope values for M0. Using first')
	if len(set(m0.header.image_defs['rescale intercept'])) > 1:
		print('Warning! Multiple rescale intercept values for M0. Using first')

	m0d = m0.get_fdata()
	pcasl4d = pcasl.get_fdata()

	# save pcasl 4D
	print('realign 4d images')
	nii = nib.as_closest_canonical(nib.Nifti1Image(pcasl4d[:,:,:,:], pcasl.affine))
	nib.save(nii,'asl4d.nii.gz')

	# realign pcasl 4d
	realigner = fsl.MCFLIRT()
	realigner.inputs.in_file='asl4d.nii.gz'
	result = realigner.run()

	#reload realigned ASL 4d
	pcaslr = nib.load(result.outputs.out_file)
	pcasl4dr = pcaslr.get_fdata()

	# take mean of ASL4d and save
	print('Registering temporal M0 mean to temporal ASL mean')
	nii = nib.as_closest_canonical(nib.Nifti1Image(m0d[:,:,:], m0.affine))
	nib.save(nii,'m0.nii.gz')

	nii = nib.as_closest_canonical(nib.Nifti1Image(pcasl4dr.mean(3), pcaslr.affine))
	nib.save(nii,'asl3d.nii.gz')

	# Register m03d to asl3d
	flt = fsl.FLIRT()
	flt.inputs.in_file = 'm0.nii.gz'
	flt.inputs.reference = 'asl3d.nii.gz'
	flt.inputs.out_file = 'm0_to_asl3d.nii.gz'
	flt.inputs.out_matrix_file = 'm0_to_asl3d.mat'
	result = flt.run()

	# reload m0 3d registered
	m0_to_asl3d = nib.load(result.outputs.out_file)
	m0r = m0_to_asl3d.get_fdata()

	# Separate asl data into ref and label
	print('Computing CBF')
	asl_label4d = pcasl4dr[:,:,:,0::2]
	asl_ref4d = pcasl4dr[:,:,:,1::2]

	# subtract and average
	asl_sub4d = asl_label4d - asl_ref4d

	# Zero voxels where either label or ref are zero
	asl_sub4d[asl_label4d==0] = 0
	asl_sub4d[asl_ref4d==0] = 0

	asl_sub3d = asl_sub4d.mean(3)

	# Loop over slices and correct for slice delay
	cbf = np.zeros(asl_sub3d.shape)
	for k in range(asl_sub3d.shape[2]):
		ePLD = PLD+(k*SliceDelay)
		cbf_num = 6000*Lambda*asl_sub3d[:,:,k]*np.exp(ePLD/T1b)
		cbf_den = 2*alpha*T1b*m0r[:,:,k]*(1-np.exp(-LabelDuration/T1b))
		cbf[:,:,k] = cbf_num/cbf_den

	# Threshold non-physical values (This is OK because GG does it!!)
	cbf = np.nan_to_num(cbf)
	cbf[cbf<0] = 0
	cbf[cbf>300] = 300

	foutput = 'cbf.nii.gz'
	cbfnii = nib.Nifti1Image(cbf,pcaslr.affine)
	nib.save(nib.as_closest_canonical(cbfnii),foutput)
	print('CBF saved as {}'.format(foutput))

if __name__ == '__main__':
	"""
	Usage: asl.py  [-p PLD] [-s SliceDelay] [-l LabelDuration] m0 pCASL
	"""

	import sys
	import getopt

	optlist, args = getopt.getopt(sys.argv[1:], 'p:s:l:')

	print(optlist,args)

	fm0 = args[0]
	fpcasl = args[1]
	PLD=1.60
	SliceDelay=0.039
	LabelDuration=1.65

	for o,a in optlist:
		if o == '-p':
			PLD = float(a)
		if o == '-s':
			SliceDelay = float(a)
		if o == '-l':
			LabelDuration = float(a)

	print(fm0,fpcasl,PLD,SliceDelay,LabelDuration)
	make_cbf(fm0,fpcasl,PLD,SliceDelay,LabelDuration)