function all_msmts = setImgParameters(file_names)
% Function to display last image of the timelapse in current folder, prompt 
% user to crop to each plate, for each plate, find the inoculum/center, and 
% also get the mask manually
% Return the plate names, crop rectangles, centers, masks as a struct array
% which will later be filled in with the polar images.

%User should have navigated to desired folder already
%Input: file names

    %make a struct with fields including "polarims" for later
    all_msmts.plate_name = {}; 
    all_msmts.crop_rect = {};
    all_msmts.colcenter = {};
    all_msmts.petrimask = {};
    all_msmts.polarims = {};
    
    user_inp = 'y'; %Keep going through plates

    while user_inp == 'y'
        
        % check if last image is the one to read in
        prompt = 'Type a number if a number other than last image should be read: ';
        to_read = input(prompt);
        if isempty(to_read)
            % Read in the last image
            I = imread(file_names{length(file_names)});
        else
            I = imread(file_names{to_read});
        end
        figure; imshow(I);
        %Ask user how many plates to do
        prompt = 'Write number of plates to do: ';
        num_plates = input(prompt); close(gcf);
        plate_names = {};
        
        % Iterate over plates
        for platenum = 1:num_plates
            chk = 0;
            while chk == 0
                try
                    disp('Done so far: '); disp(plate_names);
                    fprintf('Crop to plate number %d and double click', platenum);

                    %Crop to the plate & name
                    [im_crop, crop_rect] = imcrop(I);
                    prompt = 'Name of plate: ';
                    plate_names{platenum} = input(prompt, 's');
                    all_msmts(platenum).plate_name = plate_names(platenum);

                    % Save crop rectangle
                    all_msmts(platenum).crop_rect = crop_rect;

                    %Find the inoculum center
                    [inoc_center, thresh] = findinoc4(im_crop);
                    all_msmts(platenum).colcenter = inoc_center;

                    %Get the mask to remove rim
                    [petri1, petri_mask] = removerim_manual(im_crop);
                    all_msmts(platenum).petrimask = petri_mask;
                    chk = 1; %move to next plate
                catch
                    disp('Error occurred with this plate. Redo.');
                    resp = input('Use a different number image from timelapse instead of last im? Y/n: ', 's');
                    if strcmpi(resp, 'y')
                        useplate = input('Type the number of the image to be used.');
                        I = imread(file_names{useplate});
                        figure; imshow(I);
                        continue;
                    else
                        disp('Redoing with same image');
                        continue;
                    end
                end
            end
            
        end
        
        %Check if user would like to restart (ie messed up the plates) or finish
        %ask user to type y/n here
        prompt = 'All plates done. Redo? y/n: ';
        user_inp = input(prompt, 's');
        if strcmpi(user_inp, 'y') 
            continue
        elseif strcmpi(user_inp, 'n') 
            close;
            break
        else
            disp('Did not get y/n, ending line drawing');
            user_inp = 'n';
            close;
        end
    end

