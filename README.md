e2198_Ca_imaging
=================

This repository contains code and data for two analyses:

* Analysis of calcium imaging data and finding the cell correspondance between the calcium imaging and the EM reconstruction.
* Cell clustering based on anatomy from EM reconstruction and physiology from calcium imaging.



Data
--------------

#### EM vs Calcium
* code/cell_mapping_verified.m
	* The cell_dict_j variable in this file contains the correspondance between EM cell IDs and calciums ROI indices checked by a human expert.

#### Data for clustering and other anatomical data

* data/cell_info_clustering.mat - OBSOLETE
	* All derived data needed for the clustering task.
	Last update: 2016/4/1

* data/soma_coords_warped_20160401cleanup.mat
	Soma center coordinates transformed to the warped coordinate system.

* data/completed_051016_measurements.csv
	Soma radius and corrected segmentation volume (omnivol) for each cellid
	Currently contains info for all 396 cells

#### Ca imaging data

###### Derived data

* data/ca_dsos.20160822.mat
	* Direction selectivity and orientation selectivity computed from the parameter fits of the calcium traces.

* data/coeffs16.20160822.mat
	* Parameter fits to the calcium traces of each ROI. Specifically, coeffs16{3,2} is the double exponential fit
	 with tau, amplitudes and crossing times as parameters. #TODO: details

###### Calcium trace data

* data/roi_data.mat
	* Information pertaining to the calcium imaging in easy to access formats, including 
		calcium traces time-aligned and grouped into individual ROIs (100*DeltaF/F, before baseline detrending).
	* Running the reorganize.m script after loading this file will generate some additional variables and statistics on these data.

* data/roi_data_stimconsts.mat
	* Some addition information regarding the stimulus in the calcium imaging that are used to generate data in roi_data.mat and
	 needed for making summary plots.

###### Raw data

* data/0**_p2_g7_DSbars_200um_export.mat
	* Calcium signals processed into ROIs but not clipped (not time-aligned).

* data/sDS 8x45deg, narrow, 4.0s.QDS
	* Stimulus description script.

* data/AVG_043_p2_g7_ZoomedOut_center.tif
	* overview image

* data/AVG_4classes_043_p2_g7_ZoomedOut_center_011.rpb
	* ROIs

* rawdata/
	* These are the raw calcium imaging data.



Clustering
--------------



Ca imaging
--------------

Ca imaging data from Briggman et al. (2011).  
The stimulus is a bright bar on a dark background.
"Each of the eight stimuli used was a bar (0.2 mm wide x 1 mm long) moving in one of eight different directions at 1 mm/sec."
So the On and Off responses to the bar edges should be 1 second apart.

128 ms/frame, 31 frames per stimulus

data/roi_data.mat contains DeltaF/F for all ROIs, and some metadata


"The imaged field of view for each stimulus presentation was 100 μm × 100 μm. We acquired data from nine such fields, arranged in a 3 × 3 mosaic, for a total recorded area of 300 μm × 300 μm."


In the "data/sDS 8x45deg, narrow, 4.0s.QDS"
Note the 6th condition (direction) has a duration of 3.5s instead of 3.95s.


In the overview image (with the blood vessel ‘V’ pattern opening to the right):
0 deg (stimulus label) is left (yellow cells in Briggman et al. 2011),
90 deg is up (red cells)
180 is right (green cells)
270 is down (magenta cells)



Polar Coordinates
--------------
###### Final unified coordinates for paper:
* 0 ~ rostral
* 270 ~ ventral


###### Coordinates in code and data:

Stimuli coords (as labeled in "data/sDS 8x45deg, narrow, 4.0s.QDS"):
* 180 ~ ventral (2an Ca ~ 170deg)
* 270 ~ rostral (7o Ca ~ 250deg)

EM (Omni) y-z coord:
* 0 = +y ~ rostral
* 90 = +z ~ ventral

