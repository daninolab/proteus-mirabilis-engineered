%% medianFilter

% function file for part 3 of pre-processing image 
% before feeding into VGG-11 U-Net for ring boundary segmentation 

% implements sliding median filter of size m x n 

% inputs: 
    % aath of saved image after parts 1 & 2 of pre-processing ('intermediate images')
    % filter size m x n

function preProcImg_fin = medianFilter(preProcImg_int_path,m,n) 

    noisyImg = imread(preProcImg_int_path);
    noisyImg = im2double(noisyImg);

    fun = @(x) median(x(:));
    preProcImg_fin = nlfilter(noisyImg,'indexed',[m n],fun);
    
end

