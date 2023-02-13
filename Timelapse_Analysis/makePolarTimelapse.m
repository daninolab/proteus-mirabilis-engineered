function all_plate_msmts = makePolarTimelapse(img_folder, all_plate_msmts)
    % This function takes in a desired folder of images (user should
    % already be in directory containing it if relative path used) and
    % creates polar images + saves the basic analysis under
    % "all_plate_msmts"/in a _scanlapse2.mat file
    % images should be tifs

    % Set constants & get list of images
    date = input('Write timelapse date as mm-dd-yy: ', 's');
    exptname = input('Write name of expt (eg lrp1): ', 's');
    oldFolder = cd(img_folder);

    a = dir('*.tif'); %getlist of images in folder
    % a = dir('*.png'); %get list of images in folder
    file_names = {a.name};

    % Set constants
    interp_val = 1000;
    dpcm = 400/2.54; %~158 pix per cm

    % Get Crop Rectangle, Masks, Colony Centers for each plate
    % These will be used as inputs for flattening. Crop the image as many times
    % as there are plates. For each one, find the center. Also, manually mask
    % out the petri rim/outside space in each cropped image.
    if ~exist('all_plate_msmts', 'var')
        % Doing for the first time
        all_plate_msmts = setImgParameters(file_names);
    end
    num_plates = length(all_plate_msmts);
    num_files = length(file_names);

    % Now the main calculations and polar image creation will occur.
    try
        % Instantiate interpolants--One per plate

        chk = input('Doing interpolants. Previous interpolants will be overwritten. Continue? Y/N: ', 's');
        if strcmpi(chk, 'y') | strcmpi(chk, 'yes')

            all_plate_msmts(1).interpolants = [];
            interpolants = {}; %Cell array to store interpolants while iterating over images
            last_im = imread(file_names{length(file_names)});

            for i = 1:num_plates
                %Make interpolant:            
                crop_rect = all_plate_msmts(i).crop_rect;
                colcenter = all_plate_msmts(i).colcenter;
                petrimask = all_plate_msmts(i).petrimask;
                [F, rhoi, thetai, dist_vec] = makeInterpolant(dpcm, interp_val, last_im, crop_rect, petrimask, colcenter);

                %Display image for sanity check
                polar_img = F(thetai, rhoi); polar_img = polar_img'; 
                figure; imshow(polar_img); waitfordoubleclick; close(gcf);

                %Save interpolant, rhoi, thetai    
                all_plate_msmts(i).interpolants = F;
                all_plate_msmts(i).rhois = rhoi;
                all_plate_msmts(i).thetais = thetai;
                all_plate_msmts(i).dist_vecs = dist_vec;

                %Clear    
                disp(strcat('Done with plate: ', all_plate_msmts(i).plate_name));

            end
        else
            disp('Not redoing interpolants.');
             %continue
        end
        
        % Save this initial version of all_plate_msmts in a .mat file in 
        % corresponding folder, use a version compatible with large variables
        filename = strcat('scanlapse_msmts_2_', date, '.mat');
        save(filename, '-v7.3');

        % Make Polarim Subdirectories
        paths = strings(1, num_plates);
        for i = 1:num_plates
            % For each plate, get the plate name, add polarims to get subdirectory
            % name, make subdirectory
            plate_name = all_plate_msmts(i).plate_name{1} ;
            subdir = strcat(plate_name, '_polarims');
            mkdir(subdir);
            paths(i) = subdir;
        end

        % Initialize Storage for Avgs/Cvs/SDs
        if ~isfield(all_plate_msmts, 'avgd')
            avgds = zeros(interp_val, num_files); %Each column will be a different image
            for i = 1:length(all_plate_msmts)
                all_plate_msmts(i).avgd = avgds;
                all_plate_msmts(i).cvs = avgds;
                all_plate_msmts(i).stdevs = avgds;
            end
            clear avgds
        end

        % Iterate over ALL Images & plates & Generate Polar Images
        % NOTE: This will generate a LOT of data (eg, 308 images * 6 plates = ~16
        % gb). Make sure save location has enough storage space

        % Set wait bar & ability to cancel
        f = waitbar(0,'Starting Image 1','Name','Generating polarims...',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(f,'canceling',0);

        for im_num = 1:num_files

            % Check for clicked Cancel button
            if getappdata(f,'canceling')
                break
            end

            % Update waitbar and message
            waitbar(im_num/num_files,f,sprintf('Starting image %d of %d',im_num, num_files))

            % Prevent redoing images already done in case I had run the loop 
            % earlier and exited out before finishing:

            if isfield(all_plate_msmts, 'imfiles') && length(all_plate_msmts(num_plates).imfiles) > im_num
                %Already made polarims for all plates for the next image, skip
                if rem(im_num, 20) == 0
                    disp('Already done with image '); disp(im_num);
                end
                continue
            end

            % If not done yet:
            im_curr = imread(file_names{im_num});

            % For each image, iterate over the plates, make polarim, save, & store
            % file name in all_plate_msmts
            for plate_num = 1:length(all_plate_msmts)  
                rhoi = all_plate_msmts(plate_num).rhois;
                thetai = all_plate_msmts(plate_num).thetais;        
                crop_rect = all_plate_msmts(plate_num).crop_rect;
                mask_curr = all_plate_msmts(plate_num).petrimask;
                F = all_plate_msmts(plate_num).interpolants;
                polar_img = makePolarIm(im_curr, ...
                    crop_rect, mask_curr, F, rhoi, thetai);


                % Save each polar image as we go in the folder & store the plate name         
                plate_name = all_plate_msmts(plate_num).plate_name;
                imname = strcat(date, '_', plate_name, '_', num2str(im_num), '.tif');
                all_plate_msmts(plate_num).imfiles{im_num} = imname;
                impath = strcat(paths(plate_num), '/', imname);
                imwrite(polar_img, impath);

                % Get basic measurements 
                double_im = im2double(polar_img);
                avgd = mean(double_im, 2);
                stdevs = std(double_im, 0, 2);
                cvs = stdevs./avgd;
                %Store
                all_plate_msmts(plate_num).avgd(:, im_num) = avgd;
                all_plate_msmts(plate_num).cvs(:, im_num) = cvs;
                all_plate_msmts(plate_num).stdevs(:, im_num) = stdevs;
            end
            
            % Resave all_plate_msmts every ~20 images
            if rem(im_num, 20)==0
                save(filename, '-v7.3');
            end
        end
    
    
        delete(f);
        clear polar_img double_im avgd im_num cvs stdevs avgd plate_name imname impath
        clear im_curr crop_rect mask_curr F rhoi thetai

        disp('Done--saving analysis');

        % Save the measurements in a .mat file in corresponding folder, use a
        % version compatible with large variables
        filename = strcat('scanlapse_msmts_2_', date, '.mat');
        save(filename, '-v7.3');
       
    catch e
        fprintf(2,'The identifier was:\n%s',e.identifier);
        fprintf(2,'There was an error! The message was:\n%s',e.message);
        disp('Will return all_plate_msmts as is.');
    end
    % Finished with all steps, return to original folder
    cd(oldFolder);
end


%%%%%%%%%%%%
% Sub Functions

function [F, rhoi, thetai, dist_vec] = makeInterpolant(dpcm, interp_val, last_im, crop_rect, petrimask, colcenter)
    %Make interpolant:
    %Grop image
    crop_im = imcrop(last_im, crop_rect);

    %Mask and turn to black
    masked_plate = getMaskedPlate(crop_im, petrimask);

    %Use center & get coordinates of nonzero pixels
    X0 = colcenter(1); Y0 = colcenter(2);
    [Y, X, v]=find(masked_plate);
    X=X-X0; Y=Y-Y0;   
    %Convert coordinates to cm
    X = X/dpcm;
    Y = Y/dpcm;
    %Get the polar coordinates for the image
    [theta, rho] = cart2pol(X, Y);
    %Get the min & max polar coordinates
    rmin = min(rho); tmin = min(theta);
    rmax = max(rho); tmax = max(theta);

    %Create interpolant & remap image to the new grid
    F = scatteredInterpolant(theta, rho, v, 'nearest');
    [rhoi, thetai] = meshgrid(linspace(rmin, rmax, interp_val), linspace(tmin, tmax, interp_val));
    %Get the distance vector
    dist_vec = mean(rhoi, 1); %?????? mean(rhoi, 2);????   
end

function masked_plate = getMaskedPlate(crop_im, mask)

    im_curr = rgb2gray(im2double(crop_im));
    [row, col] = size(im_curr);
    pad_size = 300;
    pad_im = padarray(im_curr, [pad_size pad_size], 1);
    masked_im = pad_im; masked_im(mask == 0) = 0; %remove from outside of petri dish

    %Get final (unpadded) masked plate image
    masked_plate = masked_im(pad_size+1:pad_size+row, pad_size+1:pad_size+col);
end

function [polar_img, rhoi, thetai] = makePolarIm(im_curr, crop_rect, mask_curr, F, rhoi, thetai)
        % Crop to the plate
        crop_im = imcrop(im_curr, crop_rect);

        %Mask and turn to black
        masked_plate = getMaskedPlate(crop_im, mask_curr);
        
        %Find nonzero elements
        [Y, X, v]=find(masked_plate);
        
        % Pass it to the interpolant, get the polar image
        F.Values = v;

        polar_img = F(thetai, rhoi); 
        polar_img = polar_img'; 
end
