function [Nbr_sss,Now_sss] = STANCE_new_subject(study_number,new_subject_number,OKflag)
%%% Spontaneous and Task-related Activation of Neuronally Correlated Events (STANCE) %%%
% Updates global variables and creates file paths for new subject.
%
% Jason E. Hill
% STANCE_new_subject.m      updated     2 APR 2017

if ~exist('STANCE.mat','file')
    load('..\STANCE.mat');
else
    load('STANCE.mat');
end

if nargin < 3
    OKflag = false;
end

no3arginFlag = 0;
if nargin < 2
    if length(study_number) == 3
        This_sss = study_number;
        study_number = This_sss(1); 
        new_subject_number = This_sss(2);    
        no3arginFlag = 1;     
    else
        new_subject_number = Now_sss(2);
    end
end

if nargin < 1 || isempty(new_subject_number)
    study_number = Now_sss(1);
    Nbr_subjects = STANCE_how_many_subjects(study_number);
    
    new_subject_number = Nbr_subjects + 1;
    
    cd(oldDir);
end

Temp_sss(1) = study_number;
Temp_sss(2) = new_subject_number;
Temp_sss(3) = 1;

stopFlag = false;
for  level = 2:3
    
    filepath = STANCE_genpath(Temp_sss,level);
    
    if isdir(filepath) && ~stopFlag; 
        % Construct a questdlg with two options
        query = ['This folder is not empty! Remove files in ' filepath '?'];
        choice.default = 'OK';
        if OKflag
            choice = 'OK';
        else
            choice = questdlg(query, ...
                         'YesNo Menu', ...
                         'OK','No','No');
        end
        % Handle response
        switch choice
            case 'OK'
                disp(['Erasing files in ' filepath '.'])
                cmd_rmdir(filepath);
                mkdir(filepath);
            case 'No.'
                disp(['Cannot create new study in' filepath '. One already exists there.'])
                stopFlag = true;
        end
    else
        mkdir(filepath);
    end
end

if stopFlag
    % do not update anything
else
    Now_sss = Temp_sss;
    Nbr_sss(2) = Nbr_sss(2) + 1;    
    currentDir = pwd;
    cd(STANCEroot)
    save('STANCE.mat','SPMpath','subject_labels','Nbr_subject_labels','subjectsBase','subjectsT1postfix','subjectsFuzzyPostfix','filenameMNI','filenameGM','Nbr_sss','Now_sss','STANCEroot','male_labels','female_labels','M_array','MNIgmVolume');
    cd(currentDir)    
end    
    
