% Updating Timelapse Analysis 12/30/21

% Load latest
cd(timelapseDir);
allExptData = loadLatestTimelapseData();

%% Adding scan folders for new expts
for i = 12:16
    disp(allExptData(i).files);
    prompt = 'Type scanfolder location for this expt: ';
    allExptData(i).scanfolder = input(prompt, 's');
end


%% Fix timestamps for all new 2021 timelapses

% Since we got new computer can now load the exact times from the image
% names
% Then can convert that to datetime
% And create timevec as before
for i = 12:16
    img_folder = allExptData(i).scanfolder;
    time_data = getTimesFromFileNames(img_folder);
    allExptData(i).times = time_data;
    plot(time_data);
    title(i);
    waitforbuttonpress;
end

%% Fill in IPTGs

for i = 12:length(allExptData)
    % iterate over plates
    gene = lower(allExptData(i).expt);
    for j = 1:length(allExptData(i).plates)
        % get plate name
        plate = lower(allExptData(i).plates{j});
        disp(plate);
        plate = plate(1:(length(plate)-1));
        plate = erase(plate, gene);
        plate = erase(plate, 'mm');
        plate = erase(plate, 'iptg');
        plate = erase(plate, '_');
        iptg = str2double(plate);
        disp(iptg);
        allExptData(i).iptgs(j) = iptg;
    end
end

%% Adding Traj Labels to new expts 12-16

for i = 12:16
    for j = 1:4
       % retrieve lag time
       lagtime = allExptData(i).lagtime(j);
       % retrieve sw times
       swtimes = allExptData(i).swtimes{j};
       % retrieve cons times
       constimes = allExptData(i).constimes{j};
       % retrieve times vector to put labels back into numbers
       times = allExptData(i).times;
       % make labels vector
       traj_labels = ones(size(times));
       % get lag time ind
       lag_ind = find(times==lagtime);
       if isempty(lag_ind)
           [~, lag_ind] = min(abs(times-lagtime));
       end
       % make labels up until first sw phase be 0
       traj_labels(1:lag_ind) = 0;
       
       prev_phase = lagtime; % in units of h
       prev_ind = lag_ind; % in units of timepoint
       % get length of first sw phase + lag time; label sw phase 2
       
       for k = 1:length(swtimes)
           disp(k);
           % first do the sw phase
           curr_phase = swtimes(k);
           new_bound = prev_phase + curr_phase;
           new_ind = find(times==new_bound);
           if isempty(new_ind)
               [~, new_ind] = min(abs(times-new_bound));
           end
           traj_labels(prev_ind:new_ind) = 2;
           prev_phase = new_bound;
           prev_ind = new_ind;
           % now do the cons phase
           if length(constimes)>=k
               curr_phase = constimes(k);
               new_bound = prev_phase + curr_phase;
               new_ind = find(times==new_bound);
               if isempty(new_ind)
                   [~, new_ind] = min(abs(times-new_bound));
               end
               traj_labels(prev_ind:new_ind) = 1;
               prev_phase = new_bound;
               prev_ind = new_ind;
           end
           % move to next swarm phase (if any more)
%            plot(traj_labels);
%            waitforbuttonpress;
       end
%         store
        allExptData(i).traj_labels{j} = traj_labels;
    end
end

%% To do: calculate sw phase & cons phase cvs, front dens, so we can plot

% iterate over all expts
for exptnum = 12:16 %1:length(allExptData)
    disp(exptnum);
    for platenum = 1:length(allExptData(exptnum).plates)
        % for each one, load the traj labels, the front dens traj, and the front cv
        % traj
        temptrajlabels = allExptData(exptnum).traj_labels{platenum};
        tempdenstraj = allExptData(exptnum).front_dens_trajs{platenum};
        tempcvtraj = allExptData(exptnum).front_cv_trajs{platenum};

        for phaseval = 1:2
            % where swarming = 2 and consolidation = 1
            % get rid of last consolidation phase tho
            phaselabels = bwlabel(temptrajlabels==phaseval);
            total_phases = max(phaselabels);
            if phaseval==1 & total_phases>=2
%                 fprintf('Removing a phase from expt %g plate %g\n', exptnum, platenum);
                % get rid of the last consolidation phase
                phaselabels(phaselabels==total_phases) = 0;
                total_phases = total_phases-1;
            end
            dens_vals = zeros(1, total_phases);
            cv_vals = dens_vals;
            for k = 1:max(phaselabels)
                % get the curr phase
                dens_vals(k) = mean(tempdenstraj(phaselabels==k));
                cv_vals(k) = mean(tempcvtraj(phaselabels==k));
            end
%             disp(phaseval);
%             disp(dens_vals);
%             disp(cv_vals);
            % store
            if phaseval==1
                allExptData(exptnum).cons_dens_vals{platenum} = dens_vals;
                allExptData(exptnum).cons_cv_vals{platenum} = cv_vals;
            else
                allExptData(exptnum).sw_dens_vals{platenum} = dens_vals;
                allExptData(exptnum).sw_cv_vals{platenum} = cv_vals;
            end
        end
    end
end
        
     


%% To do: calculate sw phase & cons phase cvs, front dens, so we can plot

% iterate over all expts
for exptnum = 13:13  %1:length(allExptData)
    disp(exptnum);
    for platenum = 1:4 %1:length(allExptData(exptnum).plates)
        % for each one, load the traj labels, the front dens traj, and the front cv
        % traj
        temptrajlabels = allExptData(exptnum).traj_labels{platenum};

        % get the swarm speeds calculated with polyfit method
        temptraj = allExptData(exptnum).trajs{platenum};
        if iscell(temptraj)
            temptraj(cellfun(@isempty, temptraj)) = {0};
            temptraj = cell2mat(temptraj);
        end
        temptraj(temptraj==0) = 1;
        % If i is last three expts, need to get traj in cm
        if any(temptraj>4)
           tempdist = allExptData(exptnum).dists{platenum};
           temptraj = tempdist(round(temptraj));
        end
        temptraj = movmean(temptraj, 5);
        temptimes = allExptData(exptnum).times;
        p = moving_polyfit(temptimes, temptraj, 1, 5);
        
        p(p==0) = NaN;
        p(p<0.001) = NaN;
        [bin, edges] = histcounts(p(:, 1), 2);
        tempspds = p(:, 1); 
        
        %Now, get the MEDIAN swarm speeds during each swarm phase
        phaseval = 2;
        phaselabels = bwlabel(temptrajlabels==phaseval);
        total_phases = max(phaselabels);
        spd_vals = zeros(1, total_phases);
        for k = 1:max(phaselabels)
            % get the curr phase
            spd_vals(k) = nanmean(tempspds(phaselabels==k));
        end
        disp(spd_vals);
        subplot(2, 1, 1); plot(temptraj); subplot(2, 1, 2); plot(tempspds);
        waitforbuttonpress;
%         % store
%         allExptData(exptnum).sw_med_spds{platenum} = spd_vals;

    end
end




%% Mean Swarm Front Densities

% try using traj to get average front density
% figure(1);
% clf('reset');
for exptnum = 1:length(allExptData)
    fprintf('Expt %g of %g\n', exptnum, length(allExptData));
%     if length(allExptData(exptnum).local_cv_trajs) == ...
%             length(allExptData(exptnum).plates)
%         continue;
%     end
    tic;
    for j = 1:length(allExptData(exptnum).plates)
%         if length(allExptData(exptnum).local_cv_trajs) == j
%             continue;
%         end
        if size(allExptData(exptnum).local_cv_trajs{j}, 1) == 1000
            continue;
        end

        fprintf('Plate %g of %g\n', j, length(allExptData(exptnum).plates));
        tempfolder = strcat(allExptData(exptnum).exptfolder, '/', allExptData(exptnum).plates{j}, '_polarims');
    %     disp(tempfolder);
        cv_trajs_all = zeros(1000, allExptData(exptnum).num_files);
        for i = 1:allExptData(exptnum).num_files
            if rem(i, 50)==0
                disp(i);
            end
            imfile = allExptData(exptnum).imfiles{j}{i};
            
            if iscell(imfile)
                imfile = imfile{1};
                allExptData(exptnum).imfiles{j}{i} = imfile;
            end
            
            if exptnum > 1 & exptnum < 7
                % need to change out the 0's
                if strcmpi(imfile(4), '0')
                    imfile = strcat(imfile(1:3), imfile(5:end));
                end
                if strcmpi(imfile(1), '0')
                    imfile = imfile(2:end);
                end
            end
            
            if exptnum == 3 
                if j == 1
                    imfile = strrep(imfile, '_0_IPTG_umoD', '_umoD_0mM_IPTG');
                elseif j == 2
                    imfile = strrep(imfile, '_1mM_IPTG_umoD', '_umoD_1mM_IPTG');
                elseif j == 3
                    imfile = strrep(imfile, '_5mM_IPTG-umoD', '_umoD_5mM_IPTG');
                elseif j == 4
                    imfile = strrep(imfile, '_10mM_IPTG_umoD', '_umoD_10mM_IPTG');
                end
            end
                
            polarim = imread(fullfile(tempfolder, imfile));
            polarim = im2double(polarim);

            % Get location of colony front
            temptraj = allExptData(exptnum).trajs{j};
            if iscell(temptraj)
                temptraj(cellfun(@isempty, temptraj)) = {0};
                temptraj = cell2mat(temptraj);
            end
            temptraj(temptraj==0) = 1;
            % If i is last three expts, need to get traj in cm
            if exptnum>=9
               tempdist = allExptData(exptnum).dists{j};
               temptraj = tempdist(round(temptraj));
            end
            traj_loc = temptraj(i);
            if iscell(traj_loc)
                traj_loc = temptraj{i};
            end
            if isempty(traj_loc)
                traj_loc = dist_vec(1);
            end
            % convert from cm to dist
            dist_vec = allExptData(exptnum).dists{j};
            if iscell(dist_vec)
                dist_vec = dist_vec{end};
            end
            traj_ind = find(dist_vec>traj_loc, 1);

            % Trim polarim to traj location
%             polarim2 = polarim(1:traj_ind, :);
            polarim2 = polarim;
             % Get the mean local cv over the image from inoc to edge
            % Get average front CV with a window of width 10
            [~, cv_traj] = getRollingCV(polarim2, 10);
            cv_trajs_all(:, i) = cv_traj;

        end
        % Store the new measurements
        allExptData(exptnum).local_cv_trajs{j} = cv_trajs_all;
        
    end
    toc;
    
end

%% Fill in reached edge times

for i = 1:length(allExptData)
    for j = 1:length(allExptData(i).plates)
        if contains(lower(allExptData(i).plates{j}), 'flia')
            temp_val = allExptData(i).endimnum(j);
                temp_val = temptimes(temp_val);
                allExptData(i).reached_edge_time(j) = temp_val;
        end
    end
end
%%
for i = 1:length(allExptData)
    for j = 1:length(allExptData(i).plates)
        if ~contains(lower(allExptData(i).plates{j}), 'flia')
            if ~isempty(allExptData(i).endimnum)
                temp_val = allExptData(i).endimnum(j);
                if temp_val>length(temptimes)
                    temp_val = length(temptimes);
                end
                temp_val = temptimes(temp_val);
                allExptData(i).reached_edge_time(j) = temp_val;
            else
                allExptData(i).reached_edge_time(j) = 0;
            end
        end
    end
end

%% Save allExptData
% filename = file_data(1).name;
% save(filename, 'allExptData', '-v7.3');
disp('Saving updated version...');
filename = strcat('all_timelapse_data_', date, '.mat');
save(filename, 'allExptData', '-v7.3');
disp('Done.');