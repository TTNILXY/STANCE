function [SNR0_mean,SNR0_std,tSNR_std,tSNR_std] = estimate_SNR(B0,voxel_volume,FA)
% Estimate the mean and standard deviation for an EPI MRI scan.
% Author: Dr. Jason E. Hill (TTU)
% estimate SNR    updated     15 NOV 2016
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%

% Inputs:
% -------------
% voxel_volume: volume of voxel           [mm^3]
% B0: main magnet magnetic field strength [T]
% FA: flip angle                          [degrees]
%
% Outputs:
% -------------
% SNR0 = mean signal/thermal noise sigma ;i.e. a homogeneous region of a slice
% tSNR = temporal SNR variation of a single voxel

% Values are intepolated from Tables 2 & 3 in 
% "Comparison of physiological noise at 1.5 T, 3 T and 7 T and 
% optimization of fMRI acquisition parameters" by
% C. Triantafyllou, R.D. Hoge, G. Krueger, C.J. Wiggins, A. Potthast, 
% G.C. Wiggins,a and L.L. Walda, published in NeuroImage 26 (2005) 243–250.

% from Table 2:
% angles = [0,12,24,37,53,90];
% SNR0_70 = [0,30.68,  59.76,  88.43, 116.85, 184.19]; % B0 = 7.0 T
% SNR0_30 = [0,16.56, 34.72 47.40, 70.74, 99.70];      % B0 = 3.0 T
% SNR0_15 = [0, 8.48,  18.41,  31.54 , 42.82,  54.37]; % B0 = 1.5 T
% tSNR_70 = [0,28.12, 53.09, 73.74, 85.00, 95.40];     % B0 = 7.0 T
% tSNR_30 = [0,16.45, 30.28, 43.03, 56.00, 68.48];     % B0 = 3.0 T
% tSNR_15 = [0,9.43,  23.87, 30.60, 41.37, 42.38];     % B0 = 1.5 T

% from Table 3:
volumes = [3,6.75,12,27,48,75]

41.70 30.90  26.78  23.31  14.98  14.18 
77.33 51.43 46.33 37.94 24.30 22.36
131.70 F 5.07 71.80 F 3.62 1.54 F 0.14 67.98 F 1.58 55.63 F 3.43 0.70 F 0.14 42.44 F 2.79 35.87 F 1.66 0.63 F 0.1
190.29 F 5.38 77.95 F 2.86 2.23 F 0.12 87.80 F 3.19 65.62 F 1.90 0.89 F 0.09 49.98 F 4.73 42.67 F 2.71 0.61 F 0.2
346.17 F 9.66 87.23 F 4.62 3.84 F 0.25 137.51 F 3.52 73.14 F 3.28 1.59 F 0.11 74.19 F 5.41 58.53 F 4.93 0.78 F 0.2
484.21 F 6.67 92.06 F 5.44 5.16 F 0.33 199.14 F 7.52 82.09 F 3.04 2.21 F 0.14 123.95 F 5.02 73.85 F 5.28 1.35 F 0.1