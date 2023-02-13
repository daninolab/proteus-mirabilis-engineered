%% fullPreProcess
% function file for executing all pre-processing steps 
% before feeding images into VGG-11 U-Net for ring boundary segmentation 

% input: path to folder with images
    % e.g., '/Users/marianshaw/Documents/MATLAB/DIP_Spring2021/FinalProject/wt_polar_imgs'
% outputs: paths to subfolders with preprocessed images

function [preInt_sub, preFin_sub] = fullPreProcess(folderPath)

    % cd into folder
    cd(folderPath);
    
    % Overarching folder name
    [~,folderName,~] = fileparts(pwd); 
    
    % Make sub-folder to store imgs after parts 1 & 2 of pre-processing:
    % (after adaptive histogram equalization & Gaussian filtering)
    preInt_sub = strcat(folderName, '_preProc_int');
    if ~isfolder(preInt_sub)
        mkdir(preInt_sub);
    end
    
    % Make sub-folder to store imgs after part 3 of pre-processing:
    % (after Median filtering)
    preFin_sub = strcat(folderName, '_preProc_fin');
    if ~isfolder(preFin_sub)
        mkdir(preFin_sub);
    end
    
    % Get TIFF image files from folder
    imgList = dir('*.tif');
    
    % Get the names of each image file
    imgNames = {imgList.name};
   
    for imgNum = 1:length(imgList)
        
        imgName = imgNames{imgNum};
        origImg = imread(imgName);
        
        % Parts 1 & 2 of pre-processing the image 
        % (adaptive histogram equalization & Gaussian filter)
        preProcImg_int = aHist_gauss(origImg);
        % Save intermediate pre-processed images to sub-folder
        preProcImg_int_path = strcat(preInt_sub, '/', imgName);
        imwrite(preProcImg_int, preProcImg_int_path, 'tif');
        
        % Part 3 of pre-processing the image 
        % (Median fitering w/ 10x10 filter)
        preProcImg_fin = medianFilter(preProcImg_int_path,10,10);
        % Save final pre-processed images to sub-folder
        preProcImg_fin_path = strcat(preFin_sub, '/', imgName);
        imwrite(preProcImg_fin, preProcImg_fin_path, 'tif');
    end
end