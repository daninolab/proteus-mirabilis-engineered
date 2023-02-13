%% Function file for working with any file type within original folder
% output files will all be uint8 TIFF images, with 3 channels
% (written on 2-8-21, MS)

function tifNames = convertFiles(origFolder)

    % cd into folder
    cd(origFolder);

    % read files from folder
    fileList = dir();
    fileNames = {fileList.name}; % list of file names

    % first, identify non-image files to skip over
    % set up array to keep track of non-image files
    notImg = [];

    for fNum = 1:length(fileList)
        try
            isImg = imread(fileNames{fNum}); % try to read image
        catch
            warning(sprintf('%s (file # %d) is not an image.',fileNames{fNum},fNum));
            notImg = [notImg fNum];
        end
    end

    % create new list just containing image file names
    imgNames = fileNames;
    imgNames(notImg) = [];

    % 2nd, make sure images have only 3 channels and are TIFFs
    for iNum = 1:length(imgNames)
        thisImg = imread(imgNames{iNum});

        % first, check # of color channels (want only 3...4th ch = alpha channel)
        [rows cols numChannels] = size(thisImg);
        if numChannels ~= 3
            thisImg = thisImg(:,:,1:3); % extract & save only the 1st three channels 
            thisImg = im2uint8(thisImg); % convert from uint16 to uint8
        end

        % next, check image type 
        [~,name,ext] = fileparts(imgNames{iNum});
        if ~(strcmpi(ext,'.tif')) % if extension isn't tif  
            newFileName = [name, '.tif']
            imwrite(thisImg, newFileName, 'tif'); % save old jpgs/pngs etc. as new tif files 
        end
    end

    % final list of tif images
    tifList = dir('*.tif');
    tifNames = {tifList.name}; % show all tiff files 

end