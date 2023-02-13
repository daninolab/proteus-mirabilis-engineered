%% Supp Figure S3 Induction Curve Plots %%


%% Import data from excel file
% clear variables;
induction_file = '/Users/anjali/Dropbox/Anjali_Lab_Unshared/4-24-19 PM7002lacgfp_induction_phase.xlsx';

% call import script
import_gfp_data_8_11;

%% Organize Data

gfpdata = struct('IPTG', {}, 'IndTime', {}, 'GFP', {});
% PM7002lacgfplabels = PM7002lacgfplabels(:, 5:end); %get rid of leading NaNs
% PM7002lacgfptimevec = PM7002lacgfptimevec(4:(end-1)); %get rid of leading NaNs
GFParray = PM7002lacgfpinductionphaseS1(:, 3:end); %gets rid of the two wild type PM columns at start

[~, numwells] = size(GFParray); %60 wells
for i = 1:numwells
    gfpdata(i).IPTG = PM7002lacgfplabels(2, i);
    gfpdata(i).IndTime = PM7002lacgfplabels(1, i);
    gfpdata(i).GFP = GFParray(:, i)';
end

%% Calculate mean trajectories for each condition, standard error on the
% mean for each
iptgs_all = [gfpdata.IPTG];
indtimes_all = [gfpdata.IndTime];
iptgs = unique(iptgs_all);
indtimes = unique(indtimes_all);

conds =length(iptgs)*length(indtimes); %3 replicates per condition

% Define a new structure for storing these means
gfpavgs = struct('IPTG', {}, 'IndTime', {}, 'means', {}, ...
    'stdev', {}, 'sem', {});
count = 1; 
for i = 1:length(indtimes)
    indtime = indtimes(i); % the time of induction
    for j = 1:length(iptgs)
        iptg = iptgs(j);
        
        gfpavgs(count).IPTG = iptg;
        gfpavgs(count).IndTime = indtime;
        
        % Get the index where it matches both
        idxall = ([gfpdata.IPTG]==iptg) ...
            & ([gfpdata.IndTime]==indtime);
        
        % Retrieve all the data
        tempdata = cell2mat(({gfpdata(idxall).GFP})');
        
        % Calculate
        meandata = nanmean(tempdata);
        stdevs = nanstd(tempdata);
        sems = stdevs/sqrt(size(tempdata, 1));
        
        % Store
        gfpavgs(count).means = meandata;
        gfpavgs(count).stdev = stdevs;
        gfpavgs(count).sem = sems;
        
        count = count+1;
    end
    
end

clear count meandata stdevs sems tempdata idxall iptg indtime
%% Plot individual mean trajectories with standard error

figure; hold on;
t = tiledlayout(length(indtimes), 1);
for i = 1:length(gfpavgs)
    % Get vector
    meandata = gfpavgs(i).means;
    semdata = gfpavgs(i).sem;
    stdevdata = gfpavgs(i).stdev;
    % Get non-NAN locations
    inds = ~isnan(meandata);
    % Get time vector
    times = PM7002lacgfptimevec'; 
    
    ax = nexttile;
    set(ax,'FontName','Arial','FontSize',12,'TickDir','out','TickLength',...
    [0.02 0.02]);
    hold on;
    
    % plot
    errorbar(times(inds),meandata(inds),stdevdata(inds));
end
hold off;

%% Plot trajectories grouped by induction time w errorbars

% Use increasing color for more IPTG
colorvals = [0 0.33 0.66 1];

% figure; 
% t = tiledlayout(length(indtimes), 1);
for i = 1:length(indtimes)
    indtime = indtimes(i); % time of induction with IPTG
%     ax = nexttile;
    figure; ax = gca;
    set(ax,'FontName','Arial','FontSize',12,'TickDir','out','TickLength',...
    [0.02 0.02]);
    hold on;
    % Retrieve & plot reps for each iptg conc as dots
    for k = 1:length(iptgs)
        iptg = iptgs(k);
        tempcolor = colorvals(k);
        tempcolor = [tempcolor/5, tempcolor, tempcolor/5];
    
        % Get replicate trajectories and plot as dots
        repinds = find([gfpdata.IndTime]==indtime &...
            [gfpdata.IPTG] == iptg); %logical array
        for j = 1:length(repinds)
            % get the replicate
            repdata = gfpdata(repinds(j)).GFP;
%             plot(times, repdata, 'o', 'Color', tempcolor, ...
%                 'MarkerSize', 3);
                %'MarkerFaceColor', tempcolor, ... %to fill in the marker
        end     
    end
    
    
    %plot mean trajectory for each IPTG with error bars
    for k = 1:length(iptgs)
        iptg = iptgs(k);
        
        tempcolor = colorvals(k);
        tempcolor = [tempcolor/5, tempcolor, tempcolor/5];
        
        % Get mean trajectory and plot
        avgind = ([gfpavgs.IndTime]==indtime &...
            [gfpavgs.IPTG] == iptg);
        meandata = gfpavgs(avgind).means;
        semdata = gfpavgs(avgind).sem;
        stdevdata = gfpavgs(avgind).stdev;
        
        % Get the indices of non-NANs 
        inds = ~isnan(meandata);
        
        errorbar(times(inds),meandata(inds),semdata(inds),...
            'Color', tempcolor, 'MarkerSize', 0.5, ...
            'LineWidth', 1);         
    end
    
    % plot thick lines on top of error lines & for legend
    
    plotlines = gobjects(length(iptgs), 1);
    for k = 1:length(iptgs)
        iptg = iptgs(k);
        
        tempcolor = colorvals(k);
        tempcolor = [tempcolor/5, tempcolor, tempcolor/5];
        
        % Get mean trajectory and plot
        avgind = ([gfpavgs.IndTime]==indtime &...
            [gfpavgs.IPTG] == iptg);
        meandata = gfpavgs(avgind).means;
        
        % Get the indices of non-NANs 
        inds = ~isnan(meandata);
        
        % Plot a thicker line on top
        lname = strcat(num2str(iptg), ' mM IPTG');
        plotlines(k) = plot(times(inds),meandata(inds),...
            'Color', tempcolor, 'MarkerSize', 0.5, ...
            'LineWidth', 3, 'DisplayName', lname);
           
    end    
    % format figure
    xlim([0 20]); ylim([0 3000]); 
    temptitle = sprintf('Induced at %d h', indtime);
    title(temptitle,'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    xlabel('Time (h)', 'FontSize', 15);
    ylabel('Fluorescence (a.u.)', 'FontSize', 15);    
    legend(plotlines, 'Interpreter', 'none'); 
end
    

%% Plot using patch error bars
% Use increasing color for more IPTG
colorvals = [0 0.33 0.66 1];

% Define values for patches
error_alpha = .4;
error_width_factor = .003;


for i = 1:length(indtimes)

    indtime = indtimes(i); % time of induction with IPTG
    figure(i); ax = gca;
    set(ax,'FontName','Arial','FontSize',12,'TickDir','out','TickLength',...
    [0.02 0.02]);
    hold on;
 
    % store the errors for each line
    errorvals = [];
    
    %plot mean trajectory for each IPTG as line
    for k = 1:length(iptgs)
        iptg = iptgs(k);
        
        tempcolor = [colorvals(k)/5, colorvals(k), colorvals(k)/5];

        % Get mean trajectory and plot
        avgind = ([gfpavgs.IndTime]==indtime &...
            [gfpavgs.IPTG] == iptg);
        meandata = gfpavgs(avgind).means;
        semdata = gfpavgs(avgind).sem;
        stdevdata = gfpavgs(avgind).stdev;
        
        % Get the indices of non-NANs 
        inds = ~isnan(meandata);
        
        lname = strcat(num2str(iptg), ' mM IPTG');
        plotlines(k) = plot(times(inds),meandata(inds),...
            'Color', tempcolor, ...
            'LineWidth', 3, 'DisplayName', lname); 
        
        % Store the error vals
        errorvals = [errorvals; stdevdata(inds)];
    end
    
    % format figure
    xlim([0 20]); ylim([0 3000]);
    
    % plot PATCHES
    w = diff(xlim)*error_width_factor; % for patch width
    m0 = (numel(plotlines)+1)/2; % for shift
    
    for m = 1:numel(plotlines) % iterate over each line
        for n = 1:numel(plotlines(m).XData) %for each point in the line
            % make patches for each line
            patch(plotlines(m).XData(n)+[-1 1 1 -1]*w+w*(m-m0),...  
            plotlines(m).YData(n)+[-1 -1 1 1]*errorvals(m,n), 'w', 'FaceColor', ...
            plotlines(m).Color, 'FaceAlpha', error_alpha, 'EdgeColor', 'none');
        end
    end
    
 
    temptitle = sprintf('Induced at %d h', indtime);
    title(temptitle,'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
    ax.LineWidth = 1.5;
    ax.Box = 'off';
    xlabel('Time (h)', 'FontSize', 15);
    ylabel('Fluorescence (a.u.)', 'FontSize', 15);    
    legend(plotlines, 'Interpreter', 'none'); 
end


%%
figure(2); clf('reset');
% errorbar(times(inds), meandata(inds), ...
%     stdevdata(inds), 'Color', [0.2 0.5 0.188], ...
%     'LineWidth', 3);
hold on;
testx = [0, 1, 2, 3, 4, 5, 6];
testy = [0 1 2 3 4 5 6];
testvals = [0, 0.2, 0.4, 0.6, 0.8, 1];

colorvals = [0 0.33 0.66 1];
colorvals = [0 0.4 0.6 0.8];

for i = 1:6
    tempval = testvals(i);
    tempcolor = [tempval/5, tempval, tempval/5];
    plot(testx, (testy+i), 'Color', tempcolor, 'LineWidth', 3);
    
end
%% Calculate fold changes

% For each induction hour, for each induced trajectory,
% calculate the fold change at each timepoint
% between it and uninduced, then take the peak as peak fold change, plot
% for each induced timepoint vs IPTG concentration

for i = 1:length(indtimes)

    indtime = indtimes(i); % time of induction with IPTG
    % Get the baseline trajectory for that group
    baseind = ([gfpavgs.IndTime]==indtime & ...
        [gfpavgs.IPTG] == 0);
    basetraj = gfpavgs(baseind).means;
    
    %plot mean trajectory for each IPTG as line
    for k = 1:length(iptgs)
        iptg = iptgs(k);

        % Get mean trajectory 
        avgind = ([gfpavgs.IndTime]==indtime &...
            [gfpavgs.IPTG] == iptg);
        meandata = gfpavgs(avgind).means;
        
        % Get fold change & save
        foldchng = meandata./basetraj;
        gfpavgs(avgind).normtraj = foldchng;
        % Get peak fold change & save
        pkfchange = max(foldchng);
        gfpavgs(avgind).pkfchange = pkfchange;
    end
    

end
clear baseind basetraj i k iptg indtime avgind meandata foldchng pkfchange

%% PLOT PEAK FOLD CHANGE VS IPTG CONCENTRATION
fig = figure; 
ax = gca;
set(ax,'FontName','Myriad Pro','FontSize',12,'TickDir','out','TickLength',...
[0.02 0.02]);
hold on;
for i = 1:length(indtimes)
    
    indtime = indtimes(i); % time of induction with IPTG
    
    % Get the four peak fold changes for the four iptgs
    idxs = [gfpavgs.IndTime]==indtime;
    
    xvals = [gfpavgs(idxs).IPTG];
    % In order to get the 0 IPTG value on the plot, replace that x value
    % with very small number
    xvals(1) = 1e-2;
    yvals = [gfpavgs(idxs).pkfchange];
    lname = sprintf('%d h', indtime);
    h = plot(xvals, yvals, 'LineWidth', 3, 'DisplayName', lname, 'Marker', 'o',...
        'MarkerSize', 14);
    set(h, 'MarkerFaceColor', get(h,'Color'));

end

%formatting

set(ax, 'XScale', 'log');

title('pLacGFP Induction', 'FontName', 'Myriad Pro', ...
    'FontSize', 18, 'FontWeight', 'bold', 'Interpreter', 'none');
ax.LineWidth = 1.5;
ax.Box = 'off';
xlabel('IPTG (mM)', 'FontName', 'Myriad Pro', 'FontSize', 15);
ylabel('Peak Fold Change in Fluorescence', 'FontName', 'Myriad Pro', 'FontSize', 15);    
legend('Interpreter', 'none', 'FontName', 'Myriad Pro');

%% Also show a representative plot of noise between replicates

%% Plot representative growth curve--need to import other data

%% Save the gfp data
% Save updated version
disp('Saving updated version...');
filename = strcat('gfpavgs_', date, '.mat');
save(filename, 'gfpavgs', 'times');
disp('Done.');