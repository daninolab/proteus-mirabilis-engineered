function time_data = getTimesFromFileNames(img_folder)
    % Given an input folder of raw timelapse images, retrieve the times
    % they were taken/saved/created from the filenames

    oldFolder = cd(img_folder);
    a = dir('*.tif'); %getlist of images in folder
    % a = dir('*.png'); %get list of images in folder
    file_names = {a.name};
    num_files = length(file_names);
    
    %Get time 0 from first image
    im_0 = file_names{1};
    time_0 = makeImDatevec(im_0);
    
   
    %for each image, get the elapsed time
    time_data = zeros(1, length(file_names));
    for i = 1:length(file_names)
       im_time = makeImDatevec(file_names{i});
       im_etime = etime(im_time, time_0)/3600; %returns elapsed time in h
       time_data(i) = im_etime;
    end
    %time data is a vector containing elapsed time in hours for each image
    cd(oldFolder);
end

function imdate = makeImDatevec(im_name)
    % Given the tif name, make the datevec

    startind = regexp(im_name, '2021\d*_*');
    eraseind = regexp(im_name, '_\d{3}.tif')-1;
    datestr_full = im_name(startind:eraseind);
    % we don't need the milliseconds
    datestr_full = datestr_full(1:(length(datestr_full)-2));
    datestr_full = split(datestr_full, '_');
    datestr_date = datestr_full{1};
    datestr_time = datestr_full{2};
    datestr = strcat(datestr_date, "-", datestr_time);
    imdate = datevec(datestr, 'yyyymmdd-HHMMSS');
end

