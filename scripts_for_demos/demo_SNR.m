% Jason E. Hill
% demo_SNR    updated     15 NOV 2016
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%

close all; 
clear all;
currentDir = pwd;
if strcmp(currentDir(end-16:end),'scripts_for_demos')
    cd ../
    STANCEroot = pwd;    
elseif strcmp(currentDir(end-2:end),'GUI')
    % GUI instance of initialization
    cd ../
    STANCEroot = pwd;
    cd(currentDir)
elseif strcmp(currentDir(end-5:end),'STANCE')
    STANCEroot = pwd; 
else
    hSTANCE = msgbox('Please select the STANCE directory');
    uiwait(hSTANCE);
    currPath = fileparts(mfilename('fullpath'));
    STANCEroot = uigetdir(currPath, 'Add STANCE filepath');
end 
cd(STANCEroot)
addpath(genpath(pwd));

% Load STANCE globals ...
if ~exist('STANCE.mat','file')
    STANCE_initialize_STANCE;
    load('STANCE.mat');   
else
    load('STANCE.mat'); 
end
% NOTE: Must add SPM version to filepath prior to usage
addpath(SPMpath);
if exist(spm('Dir'),'dir')
    display('o SPM installation found.')
else
    warning('SPM installation not found. Please add to MATLAB filepath or install.')
    warning('SPM8 installation: http://www.fil.ion.ucl.ac.uk/spm/software/spm8/')
    exit
end


%% Turn off ...
% ... OpenGl warning
warning('off','MATLAB:opengl:StartupBlacklistedNoSetting');
% ... finite warning
warning('off', 'MATLAB:FINITE:obsoleteFunction');
% ... NIFTI class warning when loading SPM mat files
warning('off', 'MATLAB:unknownElementsNowStruc');
warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');
warning('off', 'MATLAB:pfileOlderThanMfile');

%% Modeling SNR0 versus flip angle and VV = 14.44 mm^3

% Values are intepolated from Tables 2 & 3 in 
% "Comparison of physiological noise at 1.5 T, 3 T and 7 T and 
% optimization of fMRI acquisition parameters" by
% C. Triantafyllou, R.D. Hoge, G. Krueger, C.J. Wiggins, A. Potthast, 
% G.C. Wiggins,a and L.L. Walda, published in NeuroImage 26 (2005) 243–250.


% B0 = 7T
angles = [0,12,24,37,53,90];
SNR0_70 = [0,30.68,  59.76,  88.43, 116.85, 184.19];
aq = 0:90;
SNR0q_70 = interp1(angles,SNR0_70,aq);
% B0 = 3T
SNR0_30 = [0,16.56, 34.72 47.40, 70.74, 99.70];
SNR0q_30 = interp1(angles,SNR0_30,aq);
% B0 = 1.5T
SNR0_15 = [0, 8.48,  18.41,  31.54 , 42.82,  54.37];
SNR0q_15 = interp1(angles,SNR0_15,aq);
% plot interpolation
figure
plot(angles,SNR0_70,'o',aq,SNR0q_70,':.',angles,SNR0_30,'o',aq,SNR0q_30,':.',angles,SNR0_15,'o',aq,SNR0q_15,':.');
xlim([0 90]);
title('Interpolation of SNR_0 (V = 14.44mm^3)');
legend({'7T $\mathrm{SNR}_0$ data','7T $\mathrm{SNR}_0$ linear','3T $\mathrm{SNR}_0$ data','3T $\mathrm{SNR}_0$ linear','1.5T $\mathrm{SNR}_0$ data','1.5T $\mathrm{SNR}_0$ linear'},'Location','best','Interpreter','latex')

% B0 = 7T
angles = [0,12,24,37,53,90];
SNR0_70n = [0,30.68,  59.76,  88.43, 116.85, 184.19]/184.19;
SNR0q_70n = interp1(angles,SNR0_70n,aq,'linear');
SNRmodelfun = @(b,x) (b(1)*x + b(2)*sind(x));
mdl70 = fitnlm(angles,SNR0_70n,SNRmodelfun,[1 1])
[ynew70,ynewci70] = predict(mdl70,aq');
% B0 = 3T
SNR0_30n = [0,16.56, 34.72 47.40, 70.74, 99.70]/99.70;
SNR0q_30n = interp1(angles,SNR0_30n,aq,'linear');
mdl30 = fitnlm(angles,SNR0_30n,SNRmodelfun,[1 1])
[ynew30,ynewci30] = predict(mdl30,aq');
% B0 = 1.5T
SNR0_15n = [0,8.48,  18.41,  31.54 , 42.82,  54.37]/54.37;
SNR0q_15n = interp1(angles,SNR0_15n,aq,'linear');
mdl15 = fitnlm(angles,SNR0_15n,SNRmodelfun,[1 1])
[ynew15,ynewci15] = predict(mdl15,aq');
% plot interpolation
figure
plot(angles,SNR0_70n,'o',aq,SNR0q_70n,':.',aq,ynew70,angles,SNR0_30n,'o',aq,SNR0q_30n,':.',aq,ynew30,angles,SNR0_15n,'o',aq,SNR0q_15n,':.',aq,ynew15);
xlim([0 90]);
beta1 = 1;
beta2 = 1;
x = angles;
title('Normalized SNR_0 with fit (V = 14.44mm^3)');
legend({'7T $\mathrm{SNR}_0$ data','7T $\mathrm{SNR}_0$ linear','7T $\mathrm{SNR}_0$ fit','3T $\mathrm{SNR}_0$ data','3T $\mathrm{SNR}_0$ linear','3T $\mathrm{SNR}_0$ fit','1.5T $\mathrm{SNR}_0$ data','1.5T $\mathrm{SNR}_0$ linear','1.5T $\mathrm{SNR}_0$ fit'},'Location','best','Interpreter','latex')

% model the beta parameters across B0
B0 = [1.5,3,7];
beta2 = [0.81085, 0.48239, 0.32355];
bq = [0:0.5:9.0];
beta2q = interp1(B0,beta2,bq,'spline');
beta2modelfun = @(b,x) (b(1) + b(2)*exp(-x));
mdlbeta2 = fitnlm(B0,beta2,beta2modelfun,[1 1])
[ynew,ynewci] = predict(mdlbeta2,bq');
figure, plot(B0,beta2,'o',bq,beta2q,':.',bq,ynew);
title('Fitting regression coefficient vs. B_0 (V = 14.44mm^3)');
h=legend({'$\beta_2$','spline','fit to $\beta_1 + \beta_2*\exp^{-x}$'},'Location','best','Interpreter','latex');
set(h,'Interpreter','latex')
% Resulting model fit:
% beta2 = 0.34594 + 2.1143*exp(-B0) 
% beta1 = (1 - beta2)/90 = 0.0073 - 0.0235*exp(-B0) 
% therefore:
% (0.0073 - 0.0235*exp(-B0))*FA + (0.34594 + 2.1143*exp(-B0))*sind(FA));

%% Modeling tSNR versus flip angle and V = 14.44 mm^3

% B0 = 7T
angles = [0,12,24,37,53,90];
tSNR_70 = [0,28.12, 53.09, 73.74, 85.00, 95.40];
aq = 0:90;
tSNRq_70 = interp1(angles,tSNR_70,aq,'linear');
% B0 = 3T
tSNR_30 = [0,16.45, 30.28, 43.03, 56.00, 68.48];
tSNRq_30 = interp1(angles,tSNR_30,aq,'linear');
% B0 = 1.5T
tSNR_15 = [0,9.43,  23.87, 30.60, 41.37, 42.38];
tSNRq_15 = interp1(angles,tSNR_15,aq,'linear');
% plot interpolation
figure
plot(angles,tSNR_70,'o',aq,tSNRq_70,':.',angles,tSNR_30,'o',aq,tSNRq_30,':.',angles,tSNR_15,'o',aq,tSNRq_15,':.');
xlim([0 90]);
title('Interpolation of tSNR (V = 14.44mm^3)');
legend({'7T tSNR data','7T tSNR linear','3T tSNR data','3T tSNR linear','1.5T tSNR data','1.5T tSNR linear'},'Location','best');

% B0 = 7T
angles = [0,12,24,37,53,90];
tSNR_70n = [0,28.12, 53.09, 73.74, 85.00, 95.40]/95.40;
tSNRq_70n = interp1(angles,tSNR_70n,aq,'linear');
tSNRmodelfun = @(b,x) (b(1)*x + b(2)*sind(x));
mdl70 = fitnlm(angles,tSNR_70n,tSNRmodelfun,[1 1])
[ynew70,ynewci70] = predict(mdl70,aq');
% B0 = 3T
tSNR_30n = [0,16.45, 30.28, 43.03, 56.00, 68.48]/68.48;
tSNRq_30n = interp1(angles,tSNR_30n,aq,'linear');
mdl30 = fitnlm(angles,tSNR_30n,tSNRmodelfun,[1 1])
[ynew30,ynewci30] = predict(mdl30,aq');
% B0 = 1.5T
tSNR_15n = [0,9.43,  23.87, 30.60, 41.37, 42.38]/42.38;
tSNRq_15n = interp1(angles,tSNR_15n,aq,'linear');
mdl15 = fitnlm(angles,tSNR_15n,tSNRmodelfun,[1 1])
[ynew15,ynewci15] = predict(mdl15,aq');
% plot interpolation
figure
plot(angles,tSNR_70n,'o',aq,tSNRq_70n,':.',aq,ynew70,angles,tSNR_30n,'o',aq,tSNRq_30n,':.',aq,ynew30,angles,tSNR_15n,'o',aq,tSNRq_15n,':.',aq,ynew15);
xlim([0 90]);
beta1 = 1;
beta2 = 1;
title('Normalized tSNR with fit (V = 14.44mm^3)');
legend('7T tSNR data','7T tSNR linear','7T tSNR fit','3T tSNR data','3T tSNR linear','3T tSNR fit','1.5T tSNR data','1.5T tSNR linear','1.5T tSNR fit','Location','best','Interpreter','latex');

% model consistent with simple ~sin(FA) dependancy


%% model lambda of tSNR vs SNR0

SNRq = 1:200;
tSNRq_70 = interp1(SNR0_70,tSNR_70,SNRq);
lambdaModelfun = @(b,x) (x./sqrt(1 + b^2*x.^2));
mdl70 = fitnlm(SNR0_70,tSNR_70,lambdaModelfun,0.01)
[ynew70,ynewci70] = predict(mdl70,SNRq');
%             Estimate       SE        tStat      pValue  
%  lambda:    0.0086182    0.00023954    35.979    3.1223e-07

tSNRq_30 = interp1(SNR0_30,tSNR_30,SNRq);
mdl30 = fitnlm(SNR0_30,tSNR_30,lambdaModelfun,0.01)
[ynew30,ynewci30] = predict(mdl30,SNRq');
%             Estimate       SE        tStat      pValue  
%  lambda:    0.0107      0.00027356    39.113    2.0591e-07

tSNRq_15 = interp1(SNR0_15,tSNR_15,SNRq);
mdl15 = fitnlm(SNR0_15,tSNR_15,lambdaModelfun,0.01)
[ynew15,ynewci15] = predict(mdl15,SNRq');
%             Estimate       SE        tStat      pValue  
%  lambda:    0.012333    0.0024982    4.9367    0.0043343

% plot interpolation
figure
plot(SNR0_70,tSNR_70,'o',SNRq,tSNRq_70,':.',SNRq,ynew70,SNR0_30,tSNR_30,'o',SNRq,tSNRq_30,':.',SNRq,ynew30,SNR0_15,tSNR_15,'o',SNRq,tSNRq_15,':.',SNRq,ynew15);
xlim([0 200]);
title('tSNR vs SNR_0 with \lambda fits (V = 14.44mm^3)');
legend('7T tSNR data','7T tSNR linear','7T tSNR fit','3T tSNR data','3T tSNR linear','3T tSNR fit','1.5T tSNR data','1.5T tSNR linear','1.5T tSNR fit','Location','best','Interpreter','latex')

% model the lambda parameters across TE
% (MRI protocal)
% TR = 5400 ms, three 4-mm-thick slices with a 2-mm slice gap, 60 time points,
% FOV = 240  240 mm2, matrix = 128  128, and a 
% TE of 40 ms, 30 ms, and 20 ms for 1.5 T, 3 T and 7 T, respectively.
% The bandwidth and echo spacing were, respectively 2298 Hz/pixel and
% 0.55 ms for 1.5 T, 2604 Hz/pixel and 0.43 ms for 3 T, and 2790
% Hz/pixel and 0.4 ms for the 7 T.

TE = [40, 30, 20];
TEq = 10:50;
lambda = [0.012333, 0.0107, 0.0086182];
lambda2 = lambda.^2;
lambdaTEq = interp1(TE,lambda,TEq,'spline');
lambdaTEmodelfun = @(b,x) (b(1).*(x.^2) + b(2));
mdllambdaTE = fitnlm(TE,lambda2,lambdaTEmodelfun,[1 1])
[ynewTE,ynewtTEci] = predict(mdllambdaTE,TEq');
figure, plot(TE,lambda,'o',TEq,lambdaTEq,':.',TEq,real(sqrt(ynewTE)));
xlim([0 60]);
title('\lambda vs TE fit (V = 14.44mm^3)');
legend({'$\lambda$','spline','fit to $\sqrt{\beta_1*x^2 + \beta_2}$'},'Location','best','Interpreter','latex')
%Estimated Coefficients:
%           Estimate         SE        tStat      pValue 
%    b1    6.4143e-08    7.4249e-09    8.6389    0.073365
%    b2    5.1617e-05    8.0541e-06    6.4087    0.098542
% c_1*\Delta R_2* = 2.5326e-04
%    with T_2* = 833 ms for GM, 3% activation yield \Delta T_2* = 25 ms
% c_1 = 0.0063
% c_2 = 0.0072

% model the lambda parameters across B0
lambdaq = interp1(B0,lambda,bq,'spline');
lambdaBmodelfun = @(b,x) (b(1) - b(2).*log(x));
mdllambdaB = fitnlm(B0,lambda,lambdaBmodelfun,[1 1])
[ynewB,ynewBci] = predict(mdllambdaB,bq');
figure, plot(B0,lambda,'o',bq,lambdaq,':.',bq,ynewB);
title('\lambda vs B_0 fit (V = 14.44mm^3)');
legend({'$\lambda$','spline','fit to $\beta_1 - \beta_2*\ln{x}$'},'Location','best','Interpreter','latex')
% Resulting model fit:
%          Estimate         SE        tStat      pValue  
%    b1      0.01333    3.5591e-05    374.53    0.0016998
%    b2    0.0024157    2.7143e-05    88.997     0.007153
%
% lambda = 0.01333 - 0.0024157.*log(B0)
