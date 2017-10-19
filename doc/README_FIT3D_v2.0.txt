#
# FIT3D
# 
# Version 2.0
#
# 28.09.2015
#

The new version of FIT3D comprises two different distributions, the standard
distribution in Perl, and a new distribution in Python (H.Ibarra) that it is 
distributed in a shared risk mode.

Most of the functionalities and descriptions of how to work with FIT3D included
in the README_FIT3D.txt file for version 1.0 are valid for version 2.0, 
since we tried to follow the same scheme for input and output of files. 
However, the internal algorithms have changes sustantially. For a detailed 
description of these new algorithms we refer the user to the following articles:

Sanchez et al. 2015,  RMxA&A, in print (astro-ph)

All the Perl and Python algorithms share the same inputs and output formats,
and the same name (appart from the subfix .pl instead of .py).

#
# Installation
#

For the Perl algorithms you need to install the same libraries as
the one requested by version 1.0. Then you just need to unpack the
FIT3D_v2.0.tar.gz file, and following environmental variables:

(for csh shells):
setenv FIT3D_PATH YOUR_DIRECTORY/FIT3D/
setenv PATH ${PATH}:YOUR_DIRECTORY/FIT3D/scripts/

(for bash shells):
export FIT3D_PATH =YOUR_DIRECTORY/FIT3D/
export PATH=$PATH:YOUR_DIRECTORY/FIT3D/scripts/

For the Python version you just need to unpack the tarfile, and
add the following path:

(for csh shells):
setenv PATH ${PATH}:YOUR_DIRECTORY/FIT3D_py/

(for bash shells):
export PATH=$PATH:YOUR_DIRECTORY/FIT3D_py/



#
# List of Algorithms
#

1) Algorithms for fitting the stellar population assuming a constant velocity 
dispersion in Amstrong, constant along the spectral range, that it is fitted
by the algorithm. These algorithms
are similar in this sense to version 1.0, and are useful for spectra dominated
by the instrumental resolution or at short wavelength ranges:

auto_ssp_elines_rnd.pl
auto_ssp_elines_rnd_rss.pl
auto_ssp_elines_rnd_cube.pl

 The inputs are similar to the corresponding ones in version 1.0, although
the outputs have a slightly different format. 

#
# Required inputs:
#


USE: auto_ssp_elines_rnd.pl SPEC1.txt SSP_SFH.fits,SSP_KIN.fits OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] 
CONFIG_FILE:
redshift delta_redshift min_redshift max_redshift
sigma delta_sigma min_sigma max_sigma
Av delta_Av min_Av max_Av [Same range for all]
N_SYSTEMS
(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY
...
(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY
MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX
start_w_peak end_w_peak
wavelength_to_norm width_AA new_back_templates.fits


For details on the format of the input configuration files and some examples
please read the document README_FIT3D.txt. The main difference is that it
required two input SSP files:

a) SSP_SFH.fits is the template of fitting the stellar population to derive the SFH.

b) SSP_KIN.fits is the template of fitting the stellar population to derive the stellar
kinematics.

#
# Example of its use:
#

NOTE: go to the "few_examples" directory and run:

Example 1) 
auto_ssp_elines_rnd.pl  NGC5947.spec_5.txt miles_2.fits auto_ssp.NGC5947.cen.only.out mask_elines.txt auto_ssp_V500_several_Hb.config 1 -1 40 3850 6800 emission_lines.txt  0.02 0.0001 0.005 0.03  2.5 0.25 1.2 9.0 0.5 0.1 0.0 1.6

Example 2)

auto_ssp_elines_rnd_rss.pl CS.NGC5947.RSS.fits.gz ssp_lib.fits,ssp_lib.3.fits auto_ssp.CS.NGC5947.RSS.out mask_elines.txt auto_ssp_V500_several_Hb_SII.config 1 -5 35 3700,3850 7500,5000 emission_lines.txt 0.0197 0.001 0.005 0.03 3.07 0.5 1 12.5 0.15 0.15 0.0 3.0

2) Algorithms for fitting the stellar population assuming a constant velocity 
dispersion in Amstrong, that it is fixed, in addition to a velocity dispersion in km/s
that it is fitted by the algorithm. These algorithms have the same inputs and
output than the previous one. The main difference is that they require the
instrumental velocity dispersion in AA as an extra input, and the guess, step and range
of velocity dispersions are now in km/s:


auto_ssp_elines_rnd_sigma_inst.pl
auto_ssp_elines_rnd_rss_sigma_inst.pl
auto_ssp_elines_rnd_cube_sigma_inst.pl

It is important to take into account that most SSP-templates have already an intrinsic
resolution that should be subtracted to the one included as a parameter in these algorithms.

#
# Required parameters
#
USE: auto_ssp_elines_rnd_sigma_inst.pl SPEC1.txt SSP_SFH.fits,SSP_KIN.fits,sigma_inst OUTFILE MASK_LIST CONFIG_FILE PLOT [min max] [wmin wmax] [redshift_elines_to_mask] [input_redshift delta_redshift min_redshift max_redshift] [input_sigma delta_sigma min_sigma max_sigma] [input_Av delta_Av min_Av max_Av] 
CONFIG_FILE:
redshift delta_redshift min_redshift max_redshift
sigma delta_sigma min_sigma max_sigma (km/h)
Av delta_Av min_Av max_Av [Same range for all - magnitudes ]
N_SYSTEMS
(1) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY
...
(N) START_W END_W MASK_FILE CONFIG_FILE NPOLY MASK_FILE_POLY
MIN_DELTA_CHISQ MAX_NITER CUT_MEDIAN_FLUX
start_w_peak end_w_peak
wavelength_to_norm width_AA new_back_templates.fits


#
# Example of its use:
#

NOTE: go to the "few_examples" directory and run:

Example 1) 

auto_ssp_elines_rnd_sigma_inst.pl NGC5947.spec_5.txt ssp_lib.fits,ssp_lib.3.fits,2.4 auto_ssp.NGC5947.cen.only.out mask_elines.txt auto_ssp_V500_several_Hb.config 1 -1 40 3850 6800 emission_lines.txt 0.02 0.001 0.015 0.025 10 10 5 400 0.5 0.1 0.0 1.6

3) Algorithms to fit emission lines using a set of Gaussian functions.
They have the same input and output files than the ones in version 1.0
of FIT3D (README_FIT3D.txt). However, they provide with more reliable
errors for the derived parameters:

fit_elines_rnd.pl
kin_rss_elines_rnd.pl
kin_cube_elines_rnd.pl


#
# Required inputs
#

* fit_elines_rnd.pl SPECTRUM.TXT CONFIG MASK_LIST W_MIN W_MAX [REDEFINE_MAX] [N_MC=50] [N_LOOPS=15]  [PLOT] [SCALE_INI=0.15]

In this context "SPECTRUM.TXT" is a 1D spectrum, as described in before. The CONFIG_FILE
comprises the information of the model to fit, in this particular context the most
relevants are "eline" and "poly1d". Those entries are described in README_FIT3D.txt.

The new parameters required are:
[REDEFINE_MAX]  => A flag with the options 0/1 to redefine the range of fitted
		values according to the intensity level. It is recommended to be
		fixed to 1.

[N_MC=50] 	=> Number of MC-simulations. A value of 30 it is usually fine.
[N_LOOPS=15]  	=> Number of fitting loops performed in each MC-simuation. A value of 5
		it is usually fine.
[PLOT] 		=> A flag with the values 0 (no plot), 1 (plot in the screen), 2
		plot in a PS-file.
[SCALE_INI=0.15] = > Fractional of the range of parameters for each step in each loop
		 to explore. A value of 0.15 it is usually fine.


* kin_rss_elines_rnd.pl rss.fits CONFIG MASK_LIST W_MIN W_MAX OUTFILE [N_MC=50] [N_LOOPS=15] [PLOT]  [SCALE_INI=0.15] [guided]

[guided] => corresponds to the REDEFINE_MAX in the previous algorithm

* kin_cube_elines_rnd.pl CUBE.fits CONFIG MASK_LIST W_MIN W_MAX OUTFILE [N_MC=50] [N_LOOPS=5] [PLOT] [SCALE_INI=0.15] [PREFIX] [MEMO=0/1] [input_vel_map.fits,input_vel_mask_now,FIXED=0,1] [input_disp_map.fits,FIXED=0,1 | FWHM=2.354*sigma]

[PREFIX]    => Prefix of the output maps generated for the derived parameters.
[MEMO=0/1]  => A flag indicating whether the range of parameters is stricted to that
	    of the ones derive in the previous fitting within the cube loop.

[input_vel_map.fits,input_vel_mask_now,FIXED=0,1] => Guess velocity map, and mask for follow along
						  the fitting. If FIXED is set to 1, the velocity
						  if fixed. If it is set to 0 the velocity is fitted
						  long the guess value.

[input_disp_map.fits,FIXED=0,1 | FWHM=2.354*sigma] => Guess velocity dispersion map (or guess value),
			       	 		  to follow along the fitting. If FIXED is set to 1, the 
						  velocity if fixed. If it is set to 0 the velocity 
						  is fitted long the guess value.

The two options are useful when you want to fit a weak emission line following the kinematic
pattern of a stronger one (e.g, fitting Hdelta using the kinematics information provided by Halpha).



