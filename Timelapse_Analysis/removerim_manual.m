function [petri_no_rim, plate_mask] = removerim_manual(crop_im)
J = crop_im;
grayImage = rgb2gray(im2double(J));
[row, col] = size(grayImage);
pad_size = 300;
pad_im = padarray(grayImage, [pad_size pad_size], 1);
figure(1);
imshow(pad_im);
disp('draw a circle over the actual agar');
h = drawcircle('Color','r', 'FaceAlpha',0.7); 
waitfordoubleclick;
mask = createMask(h); %returns logical with everything in circle set to true)
boundaries = bwboundaries(mask);
masked_pad = pad_im; masked_pad(mask == 0) = 1; %remove from outside of petri dish
close(figure(1));

%reverse padding--should this be to pad_size+1: (pad_size+row)?)
petri_no_rim = masked_pad(pad_size+1:pad_size+row, pad_size+1:pad_size+col);
figure(2);
imshow(petri_no_rim); waitfordoubleclick; disp('double click image to close'); close(figure(2));

plate_mask = mask;

end