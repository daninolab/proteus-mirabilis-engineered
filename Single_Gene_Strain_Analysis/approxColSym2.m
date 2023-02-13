function [CVtraj, outIm] = approxColSym2(polarim)
    % 2-20-21 method
    
    contrast_im2 = adapthisteq(polarim);
    filtim1 = medfilt2(contrast_im2);
    
    
    h = fspecial('sobel'); 
    filtim3 = imfilter(filtim1, h);
    diffim2 = imabsdiff(imresize(filtim3, size(filtim1)), filtim1);
    BW2 = imbinarize(diffim2);
    BW2_filt = bwareafilt(BW2, [50 1000000]);
    % get rid of 'largest object' ie all the white space
    
    BW3_filt = bwareafilt(BW2, 1);
    BW4 = BW2_filt;
    BW4(BW3_filt==1) = 0;
    [B,L] = bwboundaries(BW4,'noholes');
    cvs = {};
    for i = 1:length(B)
        ycoords = B{i}(:, 2);
        meancoords = movmean(ycoords, 10);
        stdevcoords = movstd(ycoords, 10);
        cvcoords = stdevcoords./meancoords;
        cvs{i} = mean(cvcoords);
    end
    
    outIm = BW4;
    CVtraj = cell2mat(cvs);
    
    
 

end