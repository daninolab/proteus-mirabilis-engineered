%% Get Swarm Speed from Timelapses V4 %%

% This script will load a 'scanlapse_2' file with the polar (flattened) timelapse images,
% iterate over each plate to get the colony edge, and save that as the swarm trajectory. 
% Using this data, swarm speed and timing will be calculated & plotted.

%% Initialize Data Structure %%
% clear all

load_val = input("Type 1 for loading data previously created by this script or 0 to load from individual scanner datas: ");
                    
if load_val == 1
    %load timelapse_swarming.mat
    [file_data, file_names, all_expts] = initAllExpts(true);
elseif load_val == 0
    %load from individual experiments
    [file_data, file_names, all_expts] = initAllExpts(false);
else
    disp('Need to type 1 or 0');
end



%% Get Swarm Trajectories

%Iterate over all the experiments & measure swarming
for expt_num = 15:15 %1:length(all_expts)

    fprintf('Starting expt num %d\n', expt_num);
    cd(allExptData(expt_num).exptfolder);
    %Iterate over the plates
    for plate_num = 4:4 %1:length(allExptData(expt_num).plates)
        tempimnames = allExptData(expt_num).imfiles{plate_num};
        fprintf('Starting plate %d \n', plate_num);
%         
%         %don't redo expts already measured
%         if isfield(allExptData, 'trajs') && ...
%             length(allExptData(expt_num).trajs)>=plate_num
%             continue
%         end
        plate_name = allExptData(expt_num).plates{plate_num};
        subdir = strcat(plate_name, '_polarims');
        cd(subdir);
        %For each plate:
        
        % Store the averages of the bg subtracted images
        avgs = {};
        % Store the locations of the swarm edge at each timepoint
        all_trajs = {};
  
        % Choose Image To Use for Start (after all condensation gone)
        start_val = 7;
        
        if isfield(allExptData, 'start_im') && ...
            length(allExptData(expt_num).start_im)>=plate_num
            start_val = allExptData(expt_num).start_im(plate_num);
        end
        
        
        fig = figure(1);
        clf('reset');
        while true
            disp(start_val);
            tempim = imread(tempimnames{start_val}{1});
            imshow(tempim);
            resp = input('Okay? type y/n: ', 's');
        
            if strcmp(resp, 'n')
                start_val = input('Type new start val: ');
                continue
            elseif strcmp(resp, 'y')
                %Save the start number in case needed later
                allExptData(expt_num).start_im(plate_num) = start_val;
                break
            else
                disp('Invalid, repeating');
                continue
            end
        end
        clear resp;
%         close(fig);
        allExptData(expt_num).start_im(j) = start_val;
        
        % New Version 6-6
        pks = {}; pks{start_val} = 1;
        avgs = {};
        start_val = start_val + 1;
        

        
        %%%%%%%%%%%%%%
        % Iterate over the time points, subtracting the previous image from each,
        % getting the average trajectory of the image, finding the edge of
        % the colony by looking for the peak difference, & saving.
        % Change the range within which to look + the minpeakprominence
        % until satisfied.
        
        % Set defaults
        prom = 0.18;
        rng = 300;
        
        if isfield(allExptData, 'prom') && ...
            length(allExptData(expt_num).prom)>=plate_num;
        
            prom = allExptData(expt_num).prom(plate_num);
            rng = allExptData(expt_num).rng(plate_num);
        end 
%         
        fig = figure(2);
        edit_colors = false;
        colors = jet;
        
        while true
            
            % For lrp, due to misaligned images, for those images just
            % assume trajectory remains in same place
%             if expt_num == 2
%                 strt = all_expts(2).all_plate_msmts(1).grps{1}(end);
%                 stp = all_expts(2).all_plate_msmts(1).grps{7}(1);
%                 imrng = [strt:stp];
%                 clear strt stp;
%             end
            
            for i = start_val:allExptData(expt_num).num_files
                %%%%%%%%%%%%%%%%%
                % 2nd 6-5 Version: Use differences, find closest peak to
                % previous
                

                % Get the difference between this im and previous im
                temp_im = imread(tempimnames{i}{1});
                temp_im_prev = imread(tempimnames{i-1}{1});
                
                if abs(double(max(max(temp_im)))-double(max(max(temp_im_prev))))>0 | ...
                        abs(double(min(min(temp_im)))-double(min(min(temp_im_prev))))>0
%                     % find the difference
%                     diff_val = double(max(max(temp_im)))-double(max(max(temp_im_prev)));
% %                     diff_val = double(temp_im(1, 500))-double(temp_im_prev(1, 500));
%                     temp_im_prev = temp_im_prev+diff_val;
% %                     
                    % Or, instead, try imhistmatch
                    temp_im_prev = imhistmatch(temp_im_prev, temp_im);
                end
                
%                 temp_im_prev = imhistmatch(temp_im_prev, temp_im);
                diffim = imabsdiff(temp_im, temp_im_prev);
                diffim = imadjust(diffim);
                z = mean(diffim, 2);       
                %z = z';
                avgs{i} = z;

                % Determine range to look for colony peak within
                lhlim = pks{i-1};
                rhlim = min(lhlim+rng, 900); 
                if i < 25
                    lhlim = 0;
                    rhlim = 100;
                elseif i < 50                
                    lhlim = pks{i-1};
                    rhlim = lhlim+50;
                end
                rhlim = min(rhlim, 950);
                [~, locs] = findpeaks(z(1:rhlim), 'MinPeakProminence', prom, ...
                    'MinPeakWidth', 1);

                if isempty(locs)
                   %If no peaks found, use previous peak location
                   pks{i} = pks{i-1};
                else
                   if locs(end)<pks{i-1}
                       % If only peaks before last peak found, use previous
                       pks{i} = pks{i-1};
                   else
                       %First get all that are greater than previous
                       temp_locs = locs(locs>lhlim); 
                       %Then use closest to previous
                       [val, ind] = min(temp_locs-lhlim);
                       % Try using furthest for capturing initial swarm
                       % phase?
                       if 40<pks{i-1}<100
                           [val, ind] = max(temp_locs-lhlim);
                       end
                       
                       pks{i} = temp_locs(ind);
                   end
                end
                %if we didn't get a peak, use the previous
                if isempty(pks{i})
                   pks{i} = pks{i-1};
                end
                
                % for lrp, for misaligned images, use previous loc
                if expt_num == 2 && i > (134) && i < (143)
                    pks{i} = pks{i-1};
                end

            end % Finished calculating over all images of plate
            
            for i = 1:length(pks)
                if isempty(pks{i})
                    pks{i} = 0;
                end
            end
            
            % refine traj: fill in first several images with inocedge
            % location
            inoc_edge_val = mode(cell2mat(pks((start_val+3):35)));
            last_inoc_reached = find(cell2mat(pks)==inoc_edge_val, 1, 'last');
            for ii = 1:last_inoc_reached
                if pks{ii} ~= inoc_edge_val | isempty(pks{ii})
                    pks{ii} = inoc_edge_val;
                end
            end
%             for i = 1:start_val
%                 pks{i} = inoc_edge_val;
%             end
%             if any(cell2mat(pks)==1)
%                 inds = find(cell2mat(pks)==1);
%                 for i = 1:length(inds)
%                     pks{inds(i)} = inoc_edge_val;
%                 end
%             end
            
            all_trajs = pks;
            %%%%%%%%%%%

            % Visualize and decide whether satisfied
            % Before visualizing, fill in []s with 0s
            temp_avgs = avgs;
            for i = 1:length(temp_avgs)
                if isempty(temp_avgs{i})
                    temp_avgs{i} = zeros(1000, 1);
                end
            end
            temp_avgs = cell2mat(temp_avgs);
            
            

%             figure; 
            clf('reset');
            hold on; h = imagesc(temp_avgs);
            colormap(h.Parent, colors);
            title(allExptData(expt_num).plates{plate_num});
            colorbar;
            %Allow increased visualization
            if edit_colors
                disp('Edit colors if desired'); 
                colormapeditor;
                waitfordoubleclick;
                colors = h.Parent.Colormap;
                edit_colors = false;
            end
            plot(cell2mat(pks), ':', 'Color', 'w', 'LineWidth', 3); 
            hold off;
            
            resp = input('Ok? Type y or n to change. ', 's');
        
            if strcmp(resp, 'n')
                % Change prominence + range & redo calculation
                fprintf('Prominence is %d and range is %d\n', prom, rng);
                prom = input('Type new prominence: ');
                rng = input('Type new range: ');
%                 close(gcf);
                continue
            elseif strcmp(resp, 'y')
                % Save the averages for that plate
                allExptData(expt_num).avgdiffs{plate_num} = avgs;
                allExptData(expt_num).trajs{plate_num} = all_trajs;
                %Save the prominence & ranges used in case needed later
                allExptData(expt_num).prom(plate_num) = prom;
                allExptData(expt_num).rng(plate_num) = rng;
%                 close(gcf);
                % Finish this plate and move to next
                break
            else
                disp('Invalid, repeating');
%                 close(gcf);
                continue
            end
            
        end % Satisfied with plate
        disp('moving on');
        cd ../
    end % Move to next plate  
    cd ../
end % Move to next expt

clear lhlim rhlim temp_im temp_im_prev bg_im diffim z avgs all_trajs;
clear expt_num plate_num;
clear start_val resp prom rng pks temp_avgs avgs val ind temp_locs

%% Fix all the distance vectors

for i = 1:length(all_expts)
    % do things differently for lrp
    if i == 2
        continue;
    end
    
    for j = 1:length(all_expts(i).all_plate_msmts)
        rhoi = all_expts(i).all_plate_msmts(j).rhois;
        %Get the distance vector
        dist_vec = mean(rhoi, 1); 
        all_expts(i).all_plate_msmts(j).dist_vecs = dist_vec;
    end
    clear rhoi dist_vec
end

%% Get Image Timestamps

for i = 6:6 %length(file_names)
    %Fill in data for that expt
   all_expts(i).files = file_names{i};
   oldFolder = cd(file_data(i).folder);
   current = file_names{i};
   time_data = getTimes(current);
   all_expts(i).all_times = time_data;
   cd(oldFolder);
end

clear time_data;

%% Save workspace (overwrite previous versions)
tic;
save('timelapse_swarming_2.mat', '-v7.3');
toc;

%% Functions

function [file_data, file_names, all_expts] = initAllExpts(loadprevious)
    oldFolder = cd('/Volumes/Seagate Backup Plus Drive/AD/scanner_timelapses');
    if loadprevious
        %load previous & add to
        load('timelapse_swarming_2.mat');
        
        %Get all analysis files done
        file_data = dir(fullfile('**', 'scanlapse_msmts_2_*.mat'));
        file_names = {file_data.name}; %cell array of file names
        
        % Fill in Data Structure --Append new ones to all_expts
        for i = 1:length(file_names)
            %Fill in data for that expt
            if ~any(strcmp({all_expts.files}, file_names{i}))
               all_expts(i).files = file_names{i};
               cd(file_data(i).folder);
               current = file_names{i};
               temp = getData(current);
               all_expts(i).all_plate_msmts =(temp{1});
               all_expts(i).num_files = temp{2};
               cd(oldFolder);
            else
                continue
            end
        end        
        
    else
        % Initialize structure
        all_expts = struct('gene', {}, 'date', {}, 'files', {}, 'all_plate_msmts', {}, 'num_files', {}); 

        % Get the Names/Locations of Msmt Files 
        
        file_data = dir(fullfile('**', 'scanlapse_msmts_2_*.mat'));
        file_names = {file_data.name}; %cell array of file names

        % Fill in Data Structure
        for i = 1:length(file_names)
            %Fill in data for that expt
           all_expts(i).files = file_names{i};
           cd(file_data(i).folder);
           current = file_names{i};
           temp = getData(current);
           all_expts(i).all_plate_msmts =(temp{1});
           all_expts(i).num_files = temp{2};
           cd(oldFolder);
        end

    end


end

%function to extract the important variables from each file

function temp_data = getData(filename)
%dostuff
    load(filename, 'all_plate_msmts', 'num_files'); disp('loaded');
    temp_data = {};
    temp_data{1} = all_plate_msmts;
    temp_data{2} = num_files;
end

function time_data = getTimes(scanlapse_path)
    load(scanlapse_path, 'file_names', 'img_folder');
    %Get time 0 from first image
    im_0 = strcat(img_folder, '/', file_names{1});
    im_0_info = dir(im_0); 
    time_0 = datevec(im_0_info.datenum);
    
    time_data = zeros(1, length(file_names));
    
    %for each image, get the elapsed time
    for i = 1:length(file_names)
       im_path = strcat(img_folder, '/', file_names{i});
       im_info = dir(im_path); 
       im_time = datevec(im_info.datenum);
       im_etime = etime(im_time, time_0)/3600; %returns elapsed time in h
       time_data(i) = im_etime;
    end
    %time data is a vector containing elapsed time in hours for each image
end