% use polyfit to get sw speeds
% figure;
cd(timelapseDir);
allExptData = loadLatestTimelapseData();
%% test on one traj
exptnum = 1;
platenum = 1;
temptraj = allExptData(1).trajs{1};
temptraj = movmean(temptraj, 3);
temptimes = allExptData(1).times;
p = moving_polyfit(temptimes, temptraj, 1, 3);

p(p==0) = NaN;
[bin, edges] = histcounts(p(:, 1), 2);
tempspds = p(:, 1); 
tempspeed = mode(tempspds(tempspds>(edges(2))));

figure(1);
clf('reset');
tempspds2 = filloutliers(tempspds, 'linear', 'gesd');
plot(temptimes, tempspds, 'o-');
hold on;
plot(temptimes, tempspds2);
xlim([0 22]);

figure(2);
clf('reset');
plot(temptimes, (tempspds'.*temptimes)+(p(:, 2)'), 'o');
hold on;
plot(temptimes, temptraj);
xlim([0 22]);

figure(3);
clf('reset');
plot(temptimes, temptraj);
hold on;
plot(temptimes, tempspds, 'o-');
xlim([0 22]);

%% Extract Swarm Times

swtimes = {};
constimes = {};

tempspds2 = movmean(tempspds, 5);
% First get the values more than 1 std above the median
highvals = tempspds2 > (median(tempspds2)+0.5*std(tempspds2));

% Define the closest regions can be.  If they are farther away than this,
% then they will be considered as separate regions.
minSeparation = 10;
highvals = ~bwareaopen(~highvals, minSeparation);
[labeledRegions, numRegions] = bwlabel(highvals);
% get bounds
bounds = cell(1, numRegions);
for i = 1:numRegions
    bounds{i} = [find(labeledRegions==i, 1, 'first'),...
        find(labeledRegions==i, 1, 'last')];
end
figure(7);
clf('reset');
plot(tempspds2);
hold on;
for i = 1:length(bounds)
    xline(bounds{i}(1));
    xline(bounds{i}(2));
end

%%
for i = 1:length(allExptData)
    for j = 1:length(allExptData(i).plates)
        temptraj = allExptData(i).trajs{j};
        if iscell(temptraj)
            temptraj(cellfun(@isempty, temptraj)) = {0};
            temptraj = cell2mat(temptraj);
        end
        temptraj(temptraj==0) = 1;
        % If i is last three expts, need to get traj in cm
        if i>=9
           tempdist = allExptData(i).dists{j};
           temptraj = tempdist(round(temptraj));
        end
        
        temptraj = movmean(temptraj, 3);
        
        temptimes = allExptData(i).times;
        p = moving_polyfit(temptimes, temptraj, 1, 5);

        p(p==0) = NaN;
        [bin, edges] = histcounts(p(:, 1), 2);
        tempspds = p(:, 1); 
        tempspeed = mode(tempspds(tempspds>(edges(2))));

        tempspds2 = movmean(tempspds, 5);
        % First get the values more than 1 std above the median
        highvals = tempspds2 > (median(tempspds2)+0.4*std(tempspds2));

        % Define the closest regions can be.  If they are farther away than this,
        % then they will be considered as separate regions.
        minSeparation = 9;
        highvals = ~bwareaopen(~highvals, minSeparation);
        [labeledRegions, numRegions] = bwlabel(highvals);
        % get bounds
        bounds = cell(1, numRegions);
        for k = 1:numRegions
            bounds{k} = [find(labeledRegions==k, 1, 'first'),...
                find(labeledRegions==k, 1, 'last')];
        end
        
        % discard very thin regions at the beginning         
        boundsmat = cell2mat(bounds);
        boundsmat2 = reshape(cell2mat(bounds'), [length(bounds), 2]);
        
        for k = 1:size(boundsmat2, 1)
            % if row contains elements less than 10, remove
            if any(boundsmat2(k, :) < 10)
                boundsmat2(k, :) = NaN;
            end
        end
        

        % discard nanrows
        boundsmat2 = boundsmat2(all(~isnan(boundsmat2),2),:);
        
        % Discard very thin regions
        tempswtimes = diff(boundsmat2, 1, 2);
        if any(tempswtimes<5) & any(tempswtimes>10)
            % these are unlikely to have happened
            inds = tempswtimes<5;
            boundsmat2 = boundsmat2(~inds, :); % remove those rows
        end
        
        % Get a version of these bounds in units of hours
        boundsmat3 = temptimes(boundsmat2);
        
        % finally get the decided regions
        swtimes = diff(boundsmat3, 1, 2);
        % get average swarm speed
        swspds = zeros(1, length(swtimes));
        for k = 1:size(boundsmat2, 1)
            % get all of the calculated polyfit speeds for this region
            spds = tempspds(boundsmat2(k, 1):boundsmat2(k, 2));
            % Since they increase and decrease at beginning/end of swarm
            % phase, take the max
            swspds(k) = max(spds);
        end
        

        % Get the constimes--they are the difference between end of each
        % row in bounds and beginning of next row
        constimes = boundsmat3(2:end, 1) - boundsmat3(1:(end-1), 2);
        % Add in lag phase
        lagtime = boundsmat3(1);
        
%         % Store values
%         allExptData(i).swtimes{j} = swtimes;
%         allExptData(i).constimes{j} = constimes;
%         allExptData(i).swspds{j} = swspds;
%         allExptData(i).lagtime(j) = lagtime;
%         
        % Store TRAJ LABELS where 0 = lag, 1 = cons, 2 = swarm
        
        traj_phase_labels = ones(1, length(temptimes));
        % set the lag labels
        traj_phase_labels(1:boundsmat2(1)) = 0;
        for m = 1:size(boundsmat3, 1)
            % get the bounds on that swarm phase
            leftbound = boundsmat2(m, 1);
            rightbound = boundsmat2(m, 2);
            traj_phase_labels(leftbound:rightbound) = 2;
        end

        allExptData(i).traj_labels{j} = traj_phase_labels;
        
%         % Visualize Traj + Swarming on Plot
%         
%         figure(7);
%         clf('reset');
%         subplot(2, 1, 1);
%         plot(temptimes, temptraj);
%         hold on;
%         for k = 1:size(boundsmat2, 1)
%             area(temptimes(boundsmat2(k,1):boundsmat2(k,2)), ...
%                 temptraj(boundsmat2(k,1):boundsmat2(k,2)),...
%                 'FaceColor', [0.8 0.8 0.8]);
% %             xline(bounds{k}(1));
% %             xline(bounds{k}(2));
%         end
% %         plot(1:20, mean(swspds)*(1:20), 'c');
%         plot((round(lagtime):(round(lagtime)+5)), ...
%             (mean(swspds)*((round(lagtime):(round(lagtime)+5))-lagtime)), 'r');
%         titletext = sprintf(strcat(allExptData(i).plates{j},...
%             ": Mean Swarm Spd %g "), mean(swspds));
%         title(titletext, 'Interpreter', 'none');
% %         title(allExptData(i).plates{j}, 'Interpreter', 'none');
% 
%         subplot(2, 1, 2);
%         plot(temptimes, traj_phase_labels);
%         
%         waitforbuttonpress;

        
        
    end
    disp(i);
end



%% SAVE THE NEW ALLeXPTDATA

filename = strcat('all_timelapse_data_', date, '.mat');
save(filename, 'allExptData', '-v7.3');
%% Create a function to retrieve values organized by gene/by iptg

function vals = getMsmtsToPlot(allExptData, field_name, iptgs, gene_list, concat)
    vals = cell(length(gene_list), length(iptgs));

    % Iterate over allexptdata, add constimes to appropriate place in data
    % structure
    for i = 1:length(allExptData)
    % Iterate over strains; create a new fig for each? or subplot
        for j = 1:length(allExptData(i).plates)
            % Iterate over conditions; 

            % Check if this plate's iptg is in the right range
            if any(iptgs==allExptData(i).iptgs(j))
                iptgind = find(iptgs==allExptData(i).iptgs(j));
                plate_curr = lower(allExptData(i).plates{j});
                for k = 1:length(gene_list)
                    if contains(plate_curr, gene_list{k})
                        geneind = k;
                    end
                end
                % Got ind of gene to use + iptgind
                % Get relevant parameter and store in vals
%                 if concat
                    tempvals = vals(geneind, iptgind);
                    if iscell(allExptData(i).(field_name)(1))
                        tempval = allExptData(i).(field_name){j};
                        tempvals = [tempvals, tempval(:).'];
                    else
                        tempval = allExptData(i).(field_name)(j);
                        tempvals = [tempvals, tempval(:).'];
                    end
                    vals{geneind, iptgind} = tempvals;
%                 else
%                     % don't concatenate
%                     tempvals = vals{geneind, iptgind};
%                     if iscell(allExptData(i).(field_name)(1))
%                         tempval = allExptData(i).(field_name){j};
%                         tempvals = {tempvals, tempval(:).'};
%                     else
%                         tempval = allExptData(i).(field_name)(j);
%                         tempvals = {tempvals, tempval(:).'};
%                     end
%                     vals{geneind, iptgind} = tempvals;
%                 end
            end
        end
    end
end
