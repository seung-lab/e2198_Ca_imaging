# e2198_Ca_imaging
Ca imaging data from Briggman et al. (2011).  
The stimulus is a bright bar on a dark background.
"Each of the eight stimuli used was a bar (0.2 mm wide x 1 mm long) moving in one of eight different directions at 1 mm/sec."
So the On and Off responses to the bar edges should be 2 seconds apart.
128 ms/frame, 31 frames per stimulus

data/roi_data.mat contains DeltaF/F for all ROIs, and some metadata


"The imaged field of view for each stimulus presentation was 100 μm × 100 μm. We acquired data from nine such fields, arranged in a 3 × 3 mosaic, for a total recorded area of 300 μm × 300 μm."


In the "data/sDS 8x45deg, narrow, 4.0s.QDS"

Note the 6th condition (direction) has a duration of 3.5s instead of 3.95s.


In the overview image (with the blood vessel ‘V’ pattern opening to the right):

0 deg is left (yellow cells in paper),
90 deg is up (red cells)
180 is right (green cells)
270 is down (magenta cells)


In Ca:
~180 is ventral (right in Ca overview, +z in omni)
?~250 is nasal  (down in Ca overview, +y in omni)

