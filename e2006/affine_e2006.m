function [A] = affine_e2006(onsacpoints, offsacpoints)
% in n * 4 matrices. 4 = 3+1 homogeneous coord

onsac_new = onsacpoints;
onsac_new(:, 1) = 62 + 3;

offsac_new = offsacpoints;
offsac_new(:, 1) = 28 - 2;

new = [onsac_new; offsac_new];
raw = [onsacpoints; offsacpoints];

A = raw \ new;
