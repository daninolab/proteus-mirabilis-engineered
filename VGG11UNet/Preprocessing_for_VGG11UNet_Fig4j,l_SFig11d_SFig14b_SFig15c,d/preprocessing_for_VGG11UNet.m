%% preprocessing_for_VGG11UNet
% Full script for pre-processing images of P. mirabilis engineered or wt strains
% before feeding into VGG-11 U-Net model for ring boundary segmentation 

clear all;

%% For working through several folders of polar/flattened images at a time

% overarching strain folder w/ experiment folders
folderPath = '/Users/marianshaw/Dropbox/polarims_processing/cheWumoD_polarims/to_move';
cd(folderPath);

% get list & names of expt folders
exptFolders = dir('*_polarims');
exptFolderNames = {exptFolders.name};

% loop through each expt folder & do preprocessing on the polar images
for folderNum = 1:length(exptFolders)
    
    % get name of expt
    exptName = exptFolderNames{folderNum};
    exptPath = strcat(folderPath, '/', exptName);
    cd(exptPath);
    
    % perform preprocessing on all polarims (+ augmented polarims) 
    [preInt_sub, preFin_sub] = fullPreProcess(exptPath);
end



