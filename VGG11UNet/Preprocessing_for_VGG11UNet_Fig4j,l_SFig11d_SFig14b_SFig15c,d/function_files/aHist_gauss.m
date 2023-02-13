%% aHist_gauss

% function file for parts 1 & 2 of pre-processing image 
% before feeding into VGG-11 U-Net for ring boundary segmentation 

% goal: improve contrast & reduce noise

% input = image matrix from reading in image file
% output = preprocessed image

function [preProcImg_int] = aHist_gauss(origImg)
    
    % Adaptive histogram equalization TWICE to improve contrast
    aHistImg = adapthisteq(origImg);
    aHistImg = adapthisteq(aHistImg);
    
    % Low-pass Gaussian filter to reduce noise
    sigG = 2; % std dev = 2
    nG = 6*sigG + 1; % filter size nGxnG
    preProcImg_int = imgaussfilt(aHistImg,sigG,'FilterSize',nG); 
    
end