%% 02/09/23 ALL FIGURE PLOTS

% This script is for generating the main figs for the paper
% which can be generated in Matlab, formatting them
% in a unified way, and saving them.

% Clear vars if needed, need to have singleExptDir, comboExptDir, and
% timelapseDir as set variables

%%
% Set constants to use throughout:

% Save to the newrevision panels folder 
% savefigdir = '/Users/anjali/Dropbox/Anjali_Updates/_Swarming_Manuscript/_High_res_figure_copies/Individual_Panels';
savefigdir = '/Users/anjali/Dropbox/Shared_Swarming_Manuscript/NatChemBio_Reviews_Round2/Updated_Figs_Captions/New_Panels/';
% heatmap_to_use = winter;
heatmap_to_use = cbrewer('seq', 'YlGnBu', 256, 'spline');
cm_heatmap_to_use = cbrewer('div', 'RdBu', 256, 'spline');
h_colormap_to_use = cbrewer('seq', 'Blues', 256, 'spline');
% Original line colors

% linecolors = [0    0.4470    0.7410;
%     0.4660    0.6740    0.1880;
%     0.8500    0.3250    0.0980;
%     0.9290    0.6940    0.1250;
%     0.4940    0.1840    0.5560;
%     0.3010    0.7450    0.9330;
%     0.6350    0.0780    0.1840];

gene_list = {'wt', 'gfp', 'chew', 'lrp', 'flia', 'umod', 'flgm'};
% New line colors as of 08/12/21; taken from:
% https://personal.sron.nl/~pault/
% muted scheme:
% linecolors = [51, 34, 136;
%     17, 119, 51;
%     68, 170, 153;
%     170, 68, 153;
%     136, 34, 85;
%     204, 102, 119;
%     136, 204, 238];

% 'bright' color scheme
linecolors = [187, 187, 187;
    34, 136, 51;
    68, 119, 170;
    238, 102, 119;
    102, 204, 238;
    170, 51, 119;
    204, 187, 68];

linecolors = linecolors/255;
tickdirval = 'out';
ticklengthval = [0.02 0.02];
plotgrid = false;
keeplinewidths = false;
axfontsize = 24; 
smallticks = false; 
titlesize = 24;

%% Supp Figure GFP Induction Plots LOADING INFO
induction_file = '/Users/anjali/Dropbox/Anjali_Lab_Unshared/4-24-19 PM7002lacgfp_induction_phase.xlsx';

% import the data to plot
indtimes = [0, 2, 4, 6, 8];
iptgs = [0 0.1 1.0 10.0];
% load gfp avgs
load('/Users/anjali/Dropbox/Patterning_Expts_Analysis/Experiments/Single_Gene_Expts/gfpavgs_23-Aug-2021.mat');
%% Supp Figure GFP Induction Plots
% Make the first set of plots with patches
% Plot using patch error bars
% Use increasing color for more IPTG
% colorvals = [0 0.33 0.66 1];
colorvals = [0 0.4 0.6 0.8];

% Define values for patches
error_alpha = .4;
error_width_factor = .003;


for i = 1:length(indtimes)

    indtime = indtimes(i); % time of induction with IPTG
    figure(i); 
    clf('reset');
    ax = gca;
%     set(ax,'FontName','Arial','FontSize',12,'TickDir','out','TickLength',...
%     [0.02 0.02]);
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
    l = legend(plotlines, 'Interpreter', 'none', 'Location', 'northeastoutside'); 
    
    xlim([0 20]);
    xticks([0 5 10 15 20]);
    ylim([0 3000]);
    yticks([0 1000 2000 3000]);
    
    keeplinewidths = true;
    axfontsize = 24;
    smallticks = false;
    titlesize = 24;
    formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...
        keeplinewidths, axfontsize, smallticks, titlesize);

    fig = gcf;

    % Things to do for saving
    set(fig, 'WindowStyle', 'normal'); % in case it was docked
    set(fig, 'Color', 'none');
    set(fig, 'Position', [438   357   722   438]); % on laptop
    pause(0.1);
    drawnow;
    figname = sprintf('FigS3_%g_h.svg', indtime);
    export_fig(fig, fullfile(savefigdir, figname), '-nocrop', '-Painters');

    fig.Color = [1 1 1];

    
    
end

%% Supp Figure GFP  PLOT PEAK FOLD CHANGE VS IPTG CONCENTRATION

% Calculate fold changes

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



fig = figure(6);
clf('reset');
ax = gca;
set(ax,'FontName','Myriad Pro','FontSize',12,'TickDir','out','TickLength',...
[0.02 0.02]);
hold on;
markers_list = {'o', 'h', 's', 'd', '^'};
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
    h = plot(xvals, yvals, 'LineWidth', 3, 'DisplayName', lname, 'Marker', markers_list{i},...
        'MarkerSize', 14, 'Color', linecolors(2, :));
    set(h, 'MarkerFaceColor', get(h,'Color'));

end

%formatting

set(ax, 'XScale', 'log');
ax.XTick = [0.0100    0.1000    1.0000   10.0000];
ax.XTickLabel = {'0', '0.1', '1', '10'};
ax.YTick = [1 2 3 4];

title('pLacGFP Induction');
% ax.LineWidth = 1.5;
% ax.Box = 'off';
xlabel('IPTG (mM)');
ylabel('Peak Fold Change in Fluorescence');    
l = legend('Interpreter', 'none', 'Location', 'northeastoutside');
% leg_pos = l.Position;
% l.Location = 'none';
% l.Position = [leg_pos(1) leg_pos(2) leg_pos(3) 2*leg_pos(4)];
% l.Location = 'bestoutside';

l.Interpreter = 'tex';
for i = 1:length(l.String)
    l.String{i} = strcat(l.String{i}, '\newline \newline'); 
end
l.Position = [(l.Position(1)+0.01) l.Position(2) (l.Position(3) + 0.01) l.Position(4)];
l.Location = 'bestoutside';

keeplinewidths = true;
axfontsize = 24;
smallticks = false;
titlesize = 24;
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...
    keeplinewidths, axfontsize, smallticks, titlesize);

fig = gcf;

% Things to do for saving
set(fig, 'WindowStyle', 'normal'); % in case it was docked
set(fig, 'Color', 'none');
set(fig, 'Position', [438   357   722   438]); % on laptop
pause(0.1);
drawnow;
figname = 'FigS3_GFP_FoldChange.svg';
export_fig(fig, fullfile(savefigdir, figname), '-nocrop', '-Painters');

fig.Color = [1 1 1];
set(fig, 'WindowStyle', 'docked');

%% cheW Plots: Load the experiment data with right concs

savefigdir = '/Users/anjali/Dropbox/Shared_Swarming_Manuscript/NatChemBio_Reviews_Round2/Updated_Figs_Captions/New_Panels/';
% Load the chew experiment with 0.5/0.7 mM 
file_data = dir(fullfile('/Users/anjali/Dropbox/Patterning_Expts_Analysis/Experiments/Single_Gene_Expts', 'cheW_fig1_data.mat'));
if isempty(file_data)
    disp('No chew experiment data found.');
else 
    load(fullfile(file_data(1).folder, file_data(1).name), 'cheWfig1exptdata');
end

%% Figure 1d: cheW Representative Trajectories Heatmap


iptg_concs = [0 0.5 0.7 1];
chew_heatmap_trajs = cell(1, 4);

for i = 1:4
    chew_heatmap_trajs{i} = {};
end

for plate_ind = 1:23


    temptraj = cheWfig1exptdata(1).avgd{plate_ind};
    tempdist = cheWfig1exptdata(1).dists{plate_ind};
    temprad = cheWfig1exptdata(1).colrads_cm_redo(plate_ind);
    temptraj(find(tempdist>temprad, 1):end) = NaN;

    % interpolate to same scale
    xvals = 0:0.0075:4; 
    % save interpolated value for mean later
    tempvec = interp1(tempdist, temptraj, xvals);

    iptg_ind = find(iptg_concs== cheWfig1exptdata(1).iptgs(plate_ind));
    tempgrp = chew_heatmap_trajs{iptg_ind};
    tempgrp{end+1, 1} = tempvec;
    chew_heatmap_trajs{iptg_ind} = tempgrp;
end

% Get the means of each traj now
mean_trajs = zeros(4, length(xvals));
for j = 1:4
    mean_trajs(j, :) = nanmean(cell2mat(chew_heatmap_trajs{j}));
end
figure(2);
clf('reset');

imAlpha=ones(size(mean_trajs));
imAlpha(isnan(mean_trajs))=0;
imagesc(xvals, (1:1:4), mean_trajs,'AlphaData',imAlpha);
caxis([0.58 0.91]);
ax = gca;
ax.XTick = [0 1 2 3 4];
ax.YTick = [1 2 3 4];
ax.YTickLabel = {'0' '0.5' '0.7' '1'};
xlabel('Distance From Center (cm)');
ylabel('IPTG (mM)');
title('cheW Strain Mean Trajectories');
colormap(heatmap_to_use);
% let's format it it not to have any box
temp_ticklengthval = [0.02 0.02];
temp_tickdirval = 'out';
axesvis = true;
rulersvis = true;
formatHeatmapTrajPlot(ax, temp_tickdirval, temp_ticklengthval, heatmap_to_use, axesvis, rulersvis, titlesize);

% % USE THIS FOR SAVING THE FIGURE
% set(gcf, 'Color', 'none');
% set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% pause(0.1);
% % set(gcf, 'Position', [680   385   560   420]);
% set(gcf, 'Position', [680   301   453   497]);
% drawnow;
% % export_fig(ax, fullfile(savefigdir, 'Fig1_D_0pt5mM_Trajs.svg'), '-nocrop', '-Painters');
% exportgraphics(ax, fullfile(savefigdir, 'Fig1_D_Rep_trajs_heatmap.eps' ));
% set(gcf, 'Color', 'white');
% set(gcf, 'WindowStyle', 'docked');



%% Fig 1e: COMBINED Col Rad + Ring Width Plot

figure(3);
clf('reset');
yyaxis left;
% First Col radius plot:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iptgs = [0, 0.5, 0.7, 1];
% Make a cell array to store the areas at each iptg, 

% Cell array: Each row represents a different gene; each column will be a
% different IPTG concentration.
chew_rads = cell(1, 4); 
temp_iptgs = cheWfig1exptdata.iptgs;
for j = 1:4
    iptg_inds = temp_iptgs==iptgs(j);
%     tempgenes = cheWfig1exptdata(i).genes;
%     gene_inds = strcmpi(tempgenes, 'chew');
%     inds = gene_inds & iptg_inds;
    inds = iptg_inds;
    % Get all the genes and areas
    temprads = cheWfig1exptdata(1).colrads_cm_redo(inds);
    chew_rads{j} = [chew_rads{j} temprads];
end


% Create cell array of means
mean_rads = zeros(1, 4);
std_rads = mean_rads;
sems_rads = mean_rads;
for i = 1:numel(mean_rads)
    mean_rads(i) = nanmean(chew_rads{i});
    tempstd = nanstd(chew_rads{i});
    std_rads(i) = tempstd;
    sems_rads(i) = tempstd/sqrt(length(chew_rads{i}));
end


% clf('reset');
hold on;
for i = 1:4
    temprads = chew_rads{i};
    tempx = i;
    s(i) = scatter(repelem(i, length(temprads)), ...
        temprads, 150, 'filled', 'MarkerFaceColor', linecolors(3, :), ...
        'MarkerFaceAlpha', 0.7, 'XJitter', 'rand', 'XJitterWidth', 0.1);
    drawnow;
    e(i) = errorbar(i, mean_rads(i), sems_rads(i), ...
        'Marker', '_', 'MarkerSize', 50, 'CapSize', 0, 'Color', [0.3 0.3 0.3], ...
        'LineWidth', 3);
end
s(4).SizeData = 50;
s(4).SizeData = 150;

% Ring width plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis right;

chew_freqs = cell(1, 4);
chew_im_names = chew_freqs;
for i = 1:4
    chew_freqs{i} = {};
end

for i = 1:23
    pkdist = cheWfig1exptdata(1).meanpeakdist(i);
    iptg_ind = find(iptg_concs== cheWfig1exptdata(1).iptgs(i));
    tempgrp = chew_freqs{iptg_ind};
    tempgrp{end+1, 1} = pkdist;

    %Store in overarching struct
    chew_freqs{iptg_ind} = tempgrp;
    tempnames = chew_im_names{iptg_ind};
    tempnames{end+1, 1} = cheWfig1exptdata(1).imfiles{i};
    chew_im_names{iptg_ind} = tempnames;
     
end

mean_ring_vals = zeros(1, 4);
std_ring_vals = mean_ring_vals;
sems_ring_vals = mean_ring_vals;
for i = 1:numel(mean_ring_vals)
    tempvals = cell2mat(chew_freqs{i});
    mean_ring_vals(i) = nanmean(tempvals);
    tempstd = nanstd(tempvals);
    std_ring_vals(i) = tempstd;
    sems_ring_vals(i) = tempstd/sqrt(length(tempvals));
end

for i = 1:4
    tempvals = cell2mat(chew_freqs{i});
%     tempx = iptgs(i);
    tempx = i;

    s(i) = scatter(repelem(i, length(tempvals)), ...
        tempvals, 150, 'filled', 'Marker', '^', ...
        'MarkerEdgeColor', linecolors(3, :), ...
        'MarkerEdgeAlpha', 1.0, ...
        'MarkerFaceColor', linecolors(3, :), ...
        'MarkerFaceAlpha', 0, ...
        'LineWidth', 2, ...
        'XJitter', 'rand', 'XJitterWidth', 0.1);
    drawnow;
    e(i) = errorbar(i, mean_ring_vals(i), sems_ring_vals(i), ...
        'Marker', '_', 'MarkerSize', 50, 'CapSize', 0, 'Color', [0.3 0.3 0.3], ...
        'LineWidth', 3);
end

% Modify plot properties
ax = gca;
xlim([0.5 4.5]);
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'0', '0.5', '0.7', '1'};
ax.XLabel.String = 'IPTG (mM)';
title('cheW Features vs IPTG')
% set(ax, 'XScale', 'log');
% ylim([0 4.6]);



yyaxis left;
% set(ax, 'XScale', 'log');
ylim([0 8]);
title('cheW Colony Radius vs IPTG')
ax.YTick = [0 2 4 6 8];
ax.YLabel.String = 'Colony Radius (cm)';

yyaxis right;
% set(ax, 'XScale', 'log');
ylim([0 7.5]);
ax.YTick = [0 2 4 6];
ax.YLabel.String = 'Ring Width (mm)';

% Format full plot style
keeplinewidths = false;
axfontsize = 24;
smallticks = false;
titlesize = 24;
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...
    keeplinewidths, axfontsize, smallticks, titlesize);

fig = gcf;

% Things to do for saving
set(fig, 'WindowStyle', 'normal'); % in case it was docked
set(fig, 'Color', 'none');
% set(gcf, 'Position', [680   385   560   420]);
set(gcf, 'Position', [518   293   621   539]);
drawnow;
export_fig(fig, fullfile(savefigdir, 'Fig1_e_Combined.svg'), '-nocrop', '-Painters');
pause(0.1);
fig.Color = [1 1 1];
set(fig, 'WindowStyle', 'docked');

%% Figure 1f Make AUCs

numims = 23;
classbins = [0, 0.1, 0.6, 0.9, 1];
classLabels = categorical({'0', '0.5', '0.7', '1'});
% Make x and y_ord vectors

x_vec = zeros(2, numims)';
iptg_vec = zeros(1, numims)';

for imind = 1:numims
    
    x_vec(imind, 1) = cheWfig1exptdata(1).meanpeakdist(imind);
    x_vec(imind, 2) = cheWfig1exptdata(1).colrads_cm_redo(imind);
    iptg_vec(imind) = cheWfig1exptdata(1).iptgs(imind);
end

% x_vec = x_vec(:, 2:3);
y_ord = ordinal(iptg_vec, {'1', '2', '3', '4'}, [], classbins);
disp(groupcounts(y_ord));

% Fit multinomial logistic model
b_ord = mnrfit(x_vec, y_ord, 'model', 'ordinal');
pihat_ord = mnrval(b_ord,x_vec,'model','ordinal');

% Calculate multiClassAUC & store
pClass = pihat_ord;
classLabel = double(y_ord);
chewauc = multiClassAUC(pClass, classLabel);
disp(chewauc);


% Plot Confus Matreix
[~, est_cats] = max(pihat_ord, [], 2);
est_cats = ordinal(est_cats);
figure(4);
clf('reset');
% plotconfusion(y, est_cats);
cm = confusionchart(double(y_ord), double(est_cats));
cm.DiagonalColor = '#3787C0';
cm.OffDiagonalColor = '#3787C0';
% cm.OffDiagonalColor = '#08306B';
titletext = strcat("cheW Prediction Accuracy binning into ", strjoin(string(classLabels), ", "));
title({titletext, sprintf("mAUC: %g", chewauc), "cheW"});
ax = gca;
ax.FontSize = 18;


set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = strcat("1e_cheWCMs.svg");
% set(gcf, 'Position', [477   238   435   429]); % on laptop
set(gcf, 'Position', [910 333 417 461]); % on monitor
pause(0.1);
drawnow;
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
pause(0.1);
set(gcf, 'WindowStyle', 'docked');
set(gcf, 'Color', 'white');
% close(gcf);


%% Load Single Gene Expt Data
% clearvars -except singleExptDir comboExptDir timelapseDir
cd(singleExptDir);
allExptData = loadLatestSingleGeneData();



%% cheW Representative trajs at 0.5 mM IPTG
 %get the 0.5 mM replicate trajectories
reptraj_inds = [13, 18; 13, 19; 13, 20; 13, 21; 27, 5; 27, 6; 27, 7];

figure(2);
clf('reset');

s = subplot(5, 1, [2 3 4 5]);
hold on;


all_trajs = {};
xvals = 0:0.0075:4; 
hold on;
for repind = 1:7

    i = reptraj_inds(repind, 1);
    j = reptraj_inds(repind, 2);

    temptraj = allExptData(i).avgd{j};
    tempdist = allExptData(i).dists{j};
    temprad = allExptData(i).colrads_cm_redo(j);
    temptraj(find(tempdist>temprad, 1):end) = NaN;

    % save interpolated value for mean later
    tempvec = interp1(tempdist, temptraj, xvals);
    all_trajs{repind} = tempvec;
    l = plot(tempdist, temptraj, ...
                            'DisplayName', sprintf('Expt %g Im %g', i, j));
    l.Color = [0.1    0.1    0.310 0.4];

end

all_trajs = all_trajs';
all_trajs_mat = cell2mat(all_trajs);
meantraj = nanmean(all_trajs_mat);
meanline = plot(xvals, nanmean(all_trajs_mat), 'LineWidth', 3, ...
    'Color', [0    0.4470    0.7410], 'DisplayName', 'Mean Traj');
% l = legend('Location', 'best');
xlabel('Distance from Center (cm)', 'FontSize', 18);
ylabel('Intensity (a.u.)', 'FontSize', 18);


ax = gca;
xlim([0 4]);
ylim([0.5 0.9]);
ax.XTick = [0 1 2 3 4];
ax.YTick = [0.5 0.7 0.9];

keeplinewidths = true;
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
keeplinewidths = false;


subplot(5, 1, 1);

% Switching to imagesc
imAlpha=ones(size(meantraj));
imAlpha(isnan(meantraj))=0;
imagesc(xvals, 1, meantraj,'AlphaData',imAlpha);
caxis([0.58 0.91]);
ax = gca;

ax.YTickLabel = {''};
colormap(heatmap_to_use);
% let's format it it not to have any box
temp_ticklengthval = [0 0];
temp_tickdirval = 'out';
axesvis = false;
rulersvis = false;
title('cheW Strain Trajectory at 0.5mM IPTG');

formatHeatmapTrajPlot(ax, temp_tickdirval, temp_ticklengthval, heatmap_to_use, ...
    axesvis, rulersvis, titlesize);

% h = heatmap(meantraj, 'GridVisible', 'off', 'MissingDataColor', 'w', 'ColorbarVisible','off');
% % h = imagesc(meantraj);
% h.XDisplayLabels = repmat({''}, size(h.XData));
% colormap(heatmap_to_use);


% exportgraphics(h, fullfile(savefigdir, 'Fig1_D_0pt5mM_Heatmap.eps' ));

% % Export full figure
% % % Try saving this figure
% set(gcf, 'Color', 'none');
% set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% % set(gcf, 'Position', [680   385   560   420]);
% pause(0.1);
% set(gcf, 'Position', [680   301   453   497]);
% pause(0.1);
% drawnow;
% exportgraphics(gcf, fullfile(savefigdir, 'Fig1_D_Combined.eps'));









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 2 Setup: Choice of IPTGs, colors
iptgs = [0, 0.1, 1.0, 10];
gene_list = {'wt', 'gfp', 'chew', 'lrp', 'flia', 'umod', 'flgm'};


%% Fig 2: Get the Indices of the Representative Images Used in Panel 1
imnames = {'1mM IPTG pLacGFP001.tif',...
    '1mM IPTG pLaccheW002.tif',...
    'pLacfliA_IPTG_5_3002.tif',...
    'pLaclrp_5IPTG_3_004.tif',...
    '10mM IPTG pLacumoD002.tif',...
    'pLacflgM_IPTG_1_3_003.tif'};
im_genes = {'gfp', 'chew', 'flia', 'lrp', 'umod', 'flgm'};
exptinds = cell(1, 6);
iminds = cell(1, 6);
for i = 1:length(imnames)
    for j = 1:length(allExptData)
        tempimfiles = allExptData(j).imfiles;
        if any(strcmpi(imnames{i}, tempimfiles))
            exptinds{i} = [exptinds{i}, j];
            iminds{i} = [iminds{i}, find(strcmpi(imnames{i}, tempimfiles))];
        end
    end
end

% Note, the gfp image is actually from 7-24-19
exptinds{1} = 19;
iminds{1} = 17;
exptinds = cell2mat(exptinds);
iminds = cell2mat(iminds);
disp(exptinds);
disp(iminds);
%
%% Figure 2: Heatmap & Trajs For Each Rep Image
% NOTE: Modify so that this goes with the above exptinds/iminds, and remove
% this all_expts_table stuff and switch to allExptData
rulersvis = false;
for i = 1:1 %length(imnames)
%     figs(i) = figure;
    figure(i);
    clf('reset');
    
    % get indices of all trajectories for it
    imind = iminds(i);
    exptind = exptinds(i);
    
    gene_curr = allExptData(exptind).genes{imind};
    gene_ind = find(strcmpi(gene_curr, gene_list));
    title(gene_curr, 'FontSize', 18);
    
    % get the trajectory
    trajcurr = allExptData(exptind).avgd{imind};
    xvals = allExptData(exptind).dists{imind};

    % Get colony radius and use it to cut trajectory off there
    % Except for GFP & umoD (hit edge of plate)
%     if i ~= 1 && i ~= 5
        radcurr = allExptData(exptind).colrads_cm_redo(imind);
        radind = find(xvals>radcurr, 1);
        trajcurr(radind:end)=NaN;
%     else
%         radind = find(xvals>4);
%         radind = radind(1)-1;
%     end
    
    % Plot the trajectory up to 4 cm
    tempind = find(xvals>4);
    tempind = tempind(1)-1; %index of last value within colony radius
    xvals_all(i, 1:tempind) = xvals(1:tempind);
    trajs_all(i, 1:tempind) = trajcurr(1:tempind);
    
    % Plot the colormap above the line intensity plot
    ax(1) = subplot(5, 1, 1);
%     need to imagesc(x, y, traj)
    % Set the NaNs to transparent
    imAlpha=ones(size(trajcurr(1:tempind)'));
    imAlpha(isnan(trajcurr(1:tempind)'))=0;
    imagesc([0, 4], [1, 1], trajcurr(1:tempind)','AlphaData',imAlpha);
    
    caxis([0.5 1]);
    colormap(heatmap_to_use);
    title(gene_curr, 'FontSize', 18);
    
    axesvis = false;
    temp_ticklengthval = [0 0];
    temp_tickdirval = 'out';
    formatHeatmapTrajPlot(ax(1), temp_tickdirval, temp_ticklengthval, heatmap_to_use, axesvis, rulersvis, titlesize);
    
    
    ax(2) = subplot(5, 1, [2 3 4 5]);
    lines = {};
%     lines{1} = plot(xvals(1:tempind), trajcurr(1:tempind));
    lines{1} = plot(xvals(1:radind), trajcurr(1:radind), 'Color', linecolors(gene_ind, :));
    xlim([0 4]); ylim([0.5 1]);
    ax(2).XTick = [0 1 2 3 4];
    ax(2).YTick = [0.5 0.75 1];
    if i == 1
        xlabeltext = 'Distance from Center (cm)';
        ylabeltext = 'Intensity (a.u.)';
        xlabel(xlabeltext);
        ylabel(ylabeltext);
    end
    titletext = '';
    linenames = {}; linenames{1} = gene_curr;
    
%     formatSwarmPlot(ax(2), xlabeltext, ylabeltext, titletext, lines, linenames)
    keeplinewidths = false;
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
    drawnow;
    pause(0.01);
%     % Save the figure
%     set(gcf, 'Color', 'none');
%     set(gcf, 'WindowStyle', 'normal'); % in case it was docked
%     savefigname = strcat("Fig2b_Rep_Heatmap_Trajs_", gene_curr, ".svg");
%     export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
%     % in case eps is needed
%     exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
%     fprintf('Done with gene %g\n', i);
%     close(gcf);
end

%% Figure 2: FTs


imnames = {'1mM IPTG pLacGFP001.tif',...
    '1mM IPTG pLaccheW002.tif',...
    'pLacfliA_IPTG_5_3002.tif',...
    'pLaclrp_5IPTG_3_004.tif',...
    '10mM IPTG pLacumoD002.tif',...
    'pLacflgM_IPTG_1_3_003.tif'};
im_genes = {'gfp', 'chew', 'flia', 'lrp', 'umod', 'flgm'};
exptinds = cell(1, 6);
iminds = cell(1, 6);
for i = 1:length(imnames)
    for j = 1:length(allExptData)
        tempimfiles = allExptData(j).imfiles;
        if any(strcmpi(imnames{i}, tempimfiles))
            exptinds{i} = [exptinds{i}, j];
            iminds{i} = [iminds{i}, find(strcmpi(imnames{i}, tempimfiles))];
        end
    end
end
disp(exptinds);
disp(iminds);

% Note, the original gfp image is actually from 7-24-19, but at some point
% we switched to using the 7-26-19 one for plotting/showing.
exptinds{1} = 19;
iminds{1} = 17;
exptinds = cell2mat(exptinds);
iminds = cell2mat(iminds);
% % Load Polar & Cartesian Im for Each
% for i = 1:length(imnames)
%     fig2ims_cart{i} = im2double(imread(fullfile(allExptData(exptinds(i)).scanfolder,...
%         allExptData(exptinds(i)).imfiles{iminds(i)})));
%     fig2ims_pol{i} = im2double(imread(fullfile(allExptData(exptinds(i)).exptfolder,...
%         allExptData(exptinds(i)).polfiles{iminds(i)})));
% end
%     
% % Convert to F Domain and Save the Ims
% imdir = fullfile(singleExptDir, 'Fourier_Transform_Images/Fig2_Ims');
% if ~isfolder(imdir)
%     mkdir(imdir);
% end
% cd(imdir);
% for i = 1:6
%     imind = iminds(i);
%     exptind = exptinds(i);
%     % Make the Images
%     tempim = fig2ims_cart{i};
%     petrimask = allExptData(exptind).petrimasks{imind};
%     colcent = [allExptData(exptind).colcenters(2*imind-1), ...
%                 allExptData(exptind).colcenters(2*imind)];
%     newim = getCentCropIm(tempim, petrimask, colcent);
%     [mag_ifftim, phase_ifftim, real_ifftim, imag_ifftim] = getIFFT(newim);
%     polim = fig2ims_pol{i};
%     [mag_fftim, phase_fftim, real_fftim, imag_fftim] = getFFT(polim);
%     
%     % Save the Images
%     imname = allExptData(exptind).imfiles{imind};
%     imname = erase(imname, '.tif');
%     newimname = strcat(imname, '_cartim_ifft_mag.tif');
%     imwrite(imadjust(mag_ifftim), fullfile(savefigdir, newimname));
%     newimname = strcat(imname, '_polim_fft_mag.tif');
%     imwrite(imadjust(rescale(mag_fftim)), fullfile(savefigdir, newimname));
% end
% 
% cd(singleExptDir);

%% Figure 2: FTs


imnames = {'1mM IPTG pLacGFP001.tif',...
    '1mM IPTG pLaccheW002.tif',...
    'pLacfliA_IPTG_5_3002.tif',...
    'pLaclrp_5IPTG_3_004.tif',...
    '10mM IPTG pLacumoD002.tif',...
    'pLacflgM_IPTG_1_3_003.tif'};
im_genes = {'gfp', 'chew', 'flia', 'lrp', 'umod', 'flgm'};
exptinds = cell(1, 6);
iminds = cell(1, 6);
for i = 1:length(imnames)
    for j = 1:length(allExptData)
        tempimfiles = allExptData(j).imfiles;
        if any(strcmpi(imnames{i}, tempimfiles))
            exptinds{i} = [exptinds{i}, j];
            iminds{i} = [iminds{i}, find(strcmpi(imnames{i}, tempimfiles))];
        end
    end
end


% Note, the original gfp image is actually from 7-24-19, but at some point
% we switched to using the 7-26-19 one for plotting/showing.
exptinds{1} = 19;
iminds{1} = 17;
disp(exptinds);
disp(iminds);

exptinds = cell2mat(exptinds);
iminds = cell2mat(iminds);
% Load Polar & Cartesian Im for Each
for i = 1:length(imnames)
    fig2ims_cart{i} = im2double(imread(fullfile(allExptData(exptinds(i)).scanfolder,...
        allExptData(exptinds(i)).imfiles{iminds(i)})));
    fig2ims_pol{i} = im2double(imread(fullfile(allExptData(exptinds(i)).exptfolder,...
        allExptData(exptinds(i)).polfiles{iminds(i)})));
end
    
% Iterate over the images, get the centered and cropped image
for i = 1:6
    imind = iminds(i);
    exptind = exptinds(i);
    % Make the Images
    tempim = fig2ims_cart{i};
    petrimask = allExptData(exptind).petrimasks{imind};
    colcent = [allExptData(exptind).colcenters(2*imind-1), ...
                allExptData(exptind).colcenters(2*imind)];
    newim = getCentCropIm(tempim, petrimask, colcent);

    [mag_ifftim, phase_ifftim, real_ifftim, imag_ifftim] = getIFFT(newim);
    polim = fig2ims_pol{i};
    [mag_fftim, phase_fftim, real_fftim, imag_fftim] = getFFT(polim);
    disp(allExptData(exptind).imfiles{imind});
    % the axes limit of the FT will be 1/2 of the width of the original im
    % 400 dpi = 400 dots per 2.54 cm
    ft_w = size(newim, 1)/(400/2.54);
    disp(ft_w/2);
end



%% Figure 2: Colony Radius (flgM)
% Make a cell array to store the areas at each iptg, for each gene
% Will then iterate over each experiment. For each, iterate over the
% possible iptgs, for any found, store in the appropriate row of the main
% cell array depending on gene. then generate a new array with means. then
% plot that vs iptg.

% Cell array: Each row represents a different gene; each column will be a
% different IPTG concentration.
all_rads = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            temprads = allExptData(i).colrads_cm_redo(inds);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_rads{temprow, j} = [all_rads{temprow, j} temprads(k)];
            end
        end
    end

end

% Create cell array of means
mean_rads = zeros(7, 4);
sems_rads = mean_rads;
for i = 1:numel(mean_rads)
    mean_rads(i) = nanmean(all_rads{i});
    tempstd = nanstd(all_rads{i});
    sems_rads(i) = tempstd/sqrt(length(all_rads{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
            if i == 7 % lrp only
    %             plot(tempxvals, mean_cvs(i, :), 'DisplayName', gene_list{i}, ...
    %                 'LineWidth', 1); 
                errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
            end
    end
    
end
l = legend('Location', 'northeastoutside');

ax = gca;
% formatLinePlot(ax);
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
ax.YTick = [1 2 3 4];
keeplinewidths=true;
axfontsize = 24; smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);
set(ax, 'XScale', 'log');
% ylim([0 1.1]);
title('Mean Col Radius (cm)')

% % Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "Fig2c_ColRadius.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));

%% Figure 2: CV (lrp) --using 10 pix sliding window

% Cell array: Each row represents a different gene; each column will be a
% different IPTG concentration.
all_cvs = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            tempcvs = allExptData(i).slide_CV_means(inds);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_cvs{temprow, j} = [all_cvs{temprow, j} tempcvs(k)];
            end
        end
    end

end

% Create cell array of means
mean_cvs = zeros(7, 4);
sems_cvs = mean_cvs;
for i = 1:numel(mean_cvs)
    mean_cvs(i) = nanmean(all_cvs{i});
    tempstd = nanstd(all_cvs{i});
    sems_cvs(i) = tempstd/sqrt(length(all_cvs{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_cvs(i, :), sems_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
            if i == 4
    %             plot(tempxvals, mean_cvs(i, :), 'DisplayName', gene_list{i}, ...
    %                 'LineWidth', 1); 
                errorbar(tempxvals, mean_cvs(i, :), sems_cvs(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
            end
    end
    
end
l = legend('Location', 'northeastoutside');
formatLegendText(l, gene_list, linecolors)
ax = gca;
set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
ax.YRuler.Exponent = 0;
ax.YTick = [0.004 0.006 0.008 0.01 0.012];
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);

% ylim([0 1.1]);
% title({'Mean Local CV', '(Window of 10 pixels slid across plate)'})
title('Mean Local CV');

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "Fig2c_CV.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));


%% Figure 2: FT mag (fliA)
patchind1 = 8;
patchind2 = 8;
all_mags = cell(7, 4); 
for i = 1:length(allExptData)
    if isempty(allExptData(i).scanfolder)
        continue;
    end
    % Get the iptgs
    temp_iptgs = allExptData(i).iptgs;
    for j = 1:4
        tempiptg = iptgs(j);
        inds = find(temp_iptgs==tempiptg);
        % Get all the genes and mags
        tempgenes = allExptData(i).genes(inds);
        tempmags = zeros(1, length(tempgenes));
        for k = 1:length(inds)
            imind = inds(k);
            % Get the magnitude array
            tempmag = allExptData(i).ftregmags{imind};
            % Get the magnitude at the location
            tempmags(k) = tempmag(patchind1, patchind2);
        end
        % Now that we have the mags & the genes, fill in an
        % organized array of magnitudes, with 7 rows & one column
        % per iptg
        for k = 1:length(tempgenes)
            % Use the gene to decide which row to store that area in
            % Use j (index of iptg conc) to decide which column to store it
            % in; extend that column to add the new area
            temprow = find(contains(gene_list, tempgenes(k)));
            all_mags{temprow, j} = [all_mags{temprow, j} tempmags(k)];
        end
    end
end

% Create cell array of means
mean_mags = zeros(7, 4);
sem_mags = mean_mags;
for i = 1:numel(mean_mags)
    mean_mags(i) = nanmean(all_mags{i});
    sem_mags(i) = nanstd(all_mags{i})/nanmean(all_mags{i});
end
        
% Plot
figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%                     plot(tempxvals, mean_mags(i, :), '--', 'DisplayName', gene_list{i}, ...
%                         'LineWidth', 2, 'Color', linecolors(i, :));
            errorbar(tempxvals, mean_mags(i, :), sem_mags(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
            if i == 5 % flia
%                     plot(tempxvals, mean_mags(i, :), 'DisplayName', gene_list{i}, ...
%                         'LineWidth', 3, 'Color', linecolors(i, :));   
                errorbar(tempxvals, mean_mags(i, :), sem_mags(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
            end
    end
end

l = legend('Location', 'northeastoutside');

ax = gca;
set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
ax.YLim = [35 85];
ax.YTick = [35 50 70 85];
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);
% ylim([0 1.1]);
title('Mean FT Region Magnitude at Region 8,8')

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "Fig2c_FTmag.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));

%% Fig 2: Inoc Edge Intens Ranges for umoD

% Cell array: Each row represents a different gene; each column will be a
% different IPTG concentration.
all_rads = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            temprads = allExptData(i).inoc_traj_rngs_col(inds);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_rads{temprow, j} = [all_rads{temprow, j} temprads(k)];
            end
        end
    end

end

% Create cell array of means
mean_rads = zeros(7, 4);
sems_rads = mean_rads;
for i = 1:numel(mean_rads)
    mean_rads(i) = nanmean(all_rads{i});
    tempstd = nanstd(all_rads{i});
    sems_rads(i) = tempstd/sqrt(length(all_rads{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
            if i == 6 % umod only
    %             plot(tempxvals, mean_cvs(i, :), 'DisplayName', gene_list{i}, ...
    %                 'LineWidth', 1); 
                errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
            end
    end
    
end
l = legend('Location', 'northeastoutside');

title('Inoculum Edge Distinctness');
ax = gca;

set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
% ax.YTick = [0 0.05 0.1 0.15 0.2];
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);


savefigdir = '/Users/anjali/Dropbox/Anjali_Updates/_Swarming_Manuscript/_High_res_figure_copies/Individual_Panels';
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
pause(0.1);
waitforbuttonpress;
savefigname = "Fig2_umoD_InocEdgeRange.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
set(gcf, 'Color', 'white');
set(gcf, 'WindowStyle', 'docked'); % in case it was docked

%% Fig 2: FT Disk Region Weight Val Wide (flgM)
    
all_rads = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            temprads = allExptData(i).inoc_traj_rngs_col(inds);
            temprads = allExptData(i).FT_weightvals_wide(inds, 2);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_rads{temprow, j} = [all_rads{temprow, j} temprads(k)];
            end
        end
    end

end

% Create cell array of means
mean_rads = zeros(7, 4);
sems_rads = mean_rads;
for i = 1:numel(mean_rads)
    mean_rads(i) = nanmean(all_rads{i});
    tempstd = nanstd(all_rads{i});
    sems_rads(i) = tempstd/sqrt(length(all_rads{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
            if i == 7 % umod only
    %             plot(tempxvals, mean_cvs(i, :), 'DisplayName', gene_list{i}, ...
    %                 'LineWidth', 1); 
                errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
            end
    end
    
end
l = legend('Location', 'northeastoutside');

title('FT Domain Weight of Central Region Radius 40');
ax = gca;

set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
% ax.YTick = [0 0.05 0.1 0.15 0.2];
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);


savefigdir = '/Users/anjali/Dropbox/Anjali_Updates/_Swarming_Manuscript/_High_res_figure_copies/Individual_Panels';
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
pause(0.1);
waitforbuttonpress;
savefigname = "Fig2_flgM_FT_Region_Weight.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
set(gcf, 'Color', 'white');
set(gcf, 'WindowStyle', 'docked'); % in case it was docked
    
    
%% Supp: Col radius (all 7)

% Cell array: Each row represents a different gene; each column will be a
% different IPTG concentration.
all_rads = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            temprads = allExptData(i).colrads_cm_redo(inds);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_rads{temprow, j} = [all_rads{temprow, j} temprads(k)];
            end
        end
    end

end

% Create cell array of means
mean_rads = zeros(7, 4);
sems_rads = mean_rads;
for i = 1:numel(mean_rads)
    mean_rads(i) = nanmean(all_rads{i});
    tempstd = nanstd(all_rads{i});
    sems_rads(i) = tempstd/sqrt(length(all_rads{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
%             if i == 7 % lrp only
    %             plot(tempxvals, mean_cvs(i, :), 'DisplayName', gene_list{i}, ...
    %                 'LineWidth', 1); 
                errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
%             end
    end
    
end
l = legend('Location', 'northeastoutside');

title('Mean Col Radius (cm)');
ax = gca;

set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
% ax.YTick = [0 1 2 3 4];
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);
% ylim([0 1.1]);


% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "SuppFig_SingleGene_AllColRads.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
pause(0.1);
close(gcf);

%% Supp: CV (all 7)
% Cell array: Each row represents a different gene; each column will be a
% different IPTG concentration.
all_cvs = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            tempcvs = allExptData(i).slide_CV_means(inds);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_cvs{temprow, j} = [all_cvs{temprow, j} tempcvs(k)];
            end
        end
    end

end

% Create cell array of means
mean_cvs = zeros(7, 4);
sems_cvs = mean_cvs;
for i = 1:numel(mean_cvs)
    mean_cvs(i) = nanmean(all_cvs{i});
    tempstd = nanstd(all_cvs{i});
    sems_cvs(i) = tempstd/sqrt(length(all_cvs{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_cvs(i, :), sems_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
                errorbar(tempxvals, mean_cvs(i, :), sems_cvs(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
    end
    
end
l = legend('Location', 'northeastoutside');

title({'Mean Local CV', '(Window of 10 pixels slid across plate)'});
ax = gca;
set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);
% ylim([0 1.1]);
ax.YRuler.Exponent = 0;
ax.YTick = [0.002 0.006 0.01 0.014];

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "SuppFig_SingleGene_AllCVs.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
pause(0.1);
close(gcf);


%% Supp: FT mag (all 7)

patchind1 = 8;
patchind2 = 8;
all_mags = cell(7, 4); 
for i = 1:length(allExptData)
    if isempty(allExptData(i).scanfolder)
        continue;
    end
    % Get the iptgs
    temp_iptgs = allExptData(i).iptgs;
    for j = 1:4
        tempiptg = iptgs(j);
        inds = find(temp_iptgs==tempiptg);
        % Get all the genes and mags
        tempgenes = allExptData(i).genes(inds);
        tempmags = zeros(1, length(tempgenes));
        for k = 1:length(inds)
            imind = inds(k);
            % Get the magnitude array
            tempmag = allExptData(i).ftregmags{imind};
            % Get the magnitude at the location
            tempmags(k) = tempmag(patchind1, patchind2);
        end
        % Now that we have the mags & the genes, fill in an
        % organized array of magnitudes, with 7 rows & one column
        % per iptg
        for k = 1:length(tempgenes)
            % Use the gene to decide which row to store that area in
            % Use j (index of iptg conc) to decide which column to store it
            % in; extend that column to add the new area
            temprow = find(contains(gene_list, tempgenes(k)));
            all_mags{temprow, j} = [all_mags{temprow, j} tempmags(k)];
        end
    end
end

% Create cell array of means
mean_mags = zeros(7, 4);
sem_mags = mean_mags;
for i = 1:numel(mean_mags)
    mean_mags(i) = nanmean(all_mags{i});
    sem_mags(i) = nanstd(all_mags{i})/nanmean(all_mags{i});
end
        
% Plot
figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%                     plot(tempxvals, mean_mags(i, :), '--', 'DisplayName', gene_list{i}, ...
%                         'LineWidth', 2, 'Color', linecolors(i, :));
            errorbar(tempxvals, mean_mags(i, :), sem_mags(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
            errorbar(tempxvals, mean_mags(i, :), sem_mags(i, :), 'DisplayName', gene_list{i}, ...
                'LineWidth', 3, 'Color', linecolors(i, :)); 
    end
end

l = legend('Location', 'northeastoutside');

ax = gca;
title('Mean FT Region Magnitude at Region 8,8');
set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
ax.YTick = [20 40 60 80 100];

axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
% ylim([0 1.1]);
formatLegendText(l, gene_list, linecolors);

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "SuppFig_SingleGene_AllFTMags.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
pause(0.1);
close(gcf);

%% SUPP Fig 2: Inoc Edge Intens Ranges 

% Cell array: Each row represents a different gene; each column will be a
% different IPTG concentration.
all_rads = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            temprads = allExptData(i).inoc_traj_rngs_col(inds);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_rads{temprow, j} = [all_rads{temprow, j} temprads(k)];
            end
        end
    end

end

% Create cell array of means
mean_rads = zeros(7, 4);
sems_rads = mean_rads;
for i = 1:numel(mean_rads)
    mean_rads(i) = nanmean(all_rads{i});
    tempstd = nanstd(all_rads{i});
    sems_rads(i) = tempstd/sqrt(length(all_rads{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
%             if i == 6 % umod only
    %             plot(tempxvals, mean_cvs(i, :), 'DisplayName', gene_list{i}, ...
    %                 'LineWidth', 1); 
                errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
%             end
    end
    
end
l = legend('Location', 'northeastoutside');

title('Inoculum Edge Distinctness');
ax = gca;

set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
% ax.YTick = [0 0.05 0.1 0.15 0.2];
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);


savefigdir = '/Users/anjali/Dropbox/Anjali_Updates/_Swarming_Manuscript/_High_res_figure_copies/Individual_Panels';
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
pause(0.1);
waitforbuttonpress;
savefigname = "Supp_Fig2_InocEdgeRange.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
set(gcf, 'Color', 'white');
set(gcf, 'WindowStyle', 'docked'); % in case it was docked

%% SUPP Fig 2: FT Disk Region Weight Val Wide (all 7);
    
all_rads = cell(7, 4); 
for i = 1:length(allExptData)
    if ~isempty(allExptData(i).scanfolder)
        % Get the iptgs
        temp_iptgs = allExptData(i).iptgs;
        for j = 1:4
            tempiptg = iptgs(j);
            inds = temp_iptgs==tempiptg;
            % Get all the genes and areas
            tempgenes = allExptData(i).genes(inds);
            temprads = allExptData(i).inoc_traj_rngs_col(inds);
            temprads = allExptData(i).FT_weightvals_wide(inds, 2);
            for k = 1:length(tempgenes)
                % Use the gene to decide which row to store that area in
                % Use j (index of iptg conc) to decide which column to store it
                % in; extend that column to add the new area
                temprow = find(contains(gene_list, tempgenes(k)));
                all_rads{temprow, j} = [all_rads{temprow, j} temprads(k)];
            end
        end
    end

end

% Create cell array of means
mean_rads = zeros(7, 4);
sems_rads = mean_rads;
for i = 1:numel(mean_rads)
    mean_rads(i) = nanmean(all_rads{i});
    tempstd = nanstd(all_rads{i});
    sems_rads(i) = tempstd/sqrt(length(all_rads{i}));
end


figure(1);
clf('reset');
set(gcf, 'Visible', 'on');
tempxvals = [0.01 0.1 1 10];
hold on;
for i = 1:length(gene_list)
    switch i
        case {1, 2}
%             plot(tempxvals, mean_cvs(i, :), '--', 'DisplayName', gene_list{i}, ...
%                 'LineWidth', 2);
            errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), '--', 'DisplayName', gene_list{i}, ...
                'LineWidth', 2, 'Color', linecolors(i, :));
        otherwise
%             if i == 7 % umod only
    %             plot(tempxvals, mean_cvs(i, :), 'DisplayName', gene_list{i}, ...
    %                 'LineWidth', 1); 
                errorbar(tempxvals, mean_rads(i, :), sems_rads(i, :), 'DisplayName', gene_list{i}, ...
                    'LineWidth', 3, 'Color', linecolors(i, :)); 
%             end
    end
    
end
l = legend('Location', 'northeastoutside');

title('FT Domain Weight of Central Region Radius 40');
ax = gca;

set(ax, 'XScale', 'log');
ax.XLim = [0.005 10];
ax.XTick = [0.01 0.1 1 10];
ax.XTickLabel = {'0', '0.1', '1', '10'};
% ax.YTick = [0 0.05 0.1 0.15 0.2];
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);


savefigdir = '/Users/anjali/Dropbox/Anjali_Updates/_Swarming_Manuscript/_High_res_figure_copies/Individual_Panels';
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
pause(0.1);
waitforbuttonpress;
savefigname = "Supp_Fig2_FT_Region_Weight.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
set(gcf, 'Color', 'white');
set(gcf, 'WindowStyle', 'docked'); % in case it was docked



%% Supp: WT Trajs from 04-9-19
wt_exptfile = '4-9-19_radial_analysis.mat';
exptind = find(strcmpi(wt_exptfile, {allExptData.exptfile}));
wt_iptgs = unique(allExptData(exptind).iptgs);

% Start figure
figure(1);
clf('reset');
hold on;
t = tiledlayout(4, 1);
xvals = 0:0.0075:4;

axesvis = true;
rulersvis = false;

for j = 1:length(wt_iptgs)
    ax = nexttile;
    % Iterate over the iptgs in the expt
    tempiptg = wt_iptgs(j);
    tempinds = find(allExptData(exptind).iptgs == tempiptg);
    all_trajs = {};
    for k = 1:length(tempinds)
        imind = tempinds(k);
        % Get the traj & the dist vec
        temptraj = allExptData(exptind).avgd{imind};
        temprad = allExptData(exptind).colrads_pix_redo(imind);
        tempradcm = allExptData(exptind).colrads_cm_redo(imind);
        % Interpolate the vector from 0 to 4
        % Get the prior dist vector
        tempdist = allExptData(exptind).dists{imind};
        % Make the new trajectory via interpolation (from 0 to 4 cm)
        tempvec = interp1(tempdist, temptraj, xvals);
        if isnan(tempvec(1))
            tempvec(1) = tempvec(2);
        end
        % Now reset anything outside of the colony to NaN
        radind = find((xvals>tempradcm), 1);
        tempvec(radind:end) = NaN;
        all_trajs{k} = tempvec;
    end
    
    % Convert the traj cell to a matrix
    all_trajs = all_trajs';
    temptrajs = vertcat(all_trajs);
    yvals = 1:size(temptrajs, 1);
    temptrajs = cell2mat(temptrajs);
    
    % Imagesc, replace NaNs with white
    imAlpha=ones(size(temptrajs));
    imAlpha(isnan(temptrajs))=0;
    imagesc(xvals, yvals, temptrajs,'AlphaData',imAlpha);
    % Make sure the color axis represents the same thing
    colormap(heatmap_to_use);
    caxis([0.6 1]);
    % Add Y Label as appropriate
%     ylabeltext = sprintf('%g mm IPTG', tempiptg);
    ylabeltext = tempiptg;
    ax.YLabel.String = ylabeltext;
    
    % Add colorbar to one side
    if j == 2
        c = colorbar;
        c.TickLength = 0; c.Box = 'off';
        c.AxisLocation = "out";
        c.Ticks = [0.6 0.7 0.8 0.9 1];
        c.TickLabels = {'0.6', '0.7', '0.8', '0.9', '1'};
        drawnow;
        pause(0.1);
        c.Ruler.Axle.Visible = 'off';
    end
    
    % Add X Label/ticks to the bottom most group
    if j == 1
        title('WT Trajectories at Range of IPTG Concentrations');
    end
    if j ~= 4
        axesvis = true;
        rulersvis = false;
        ax.XTick = [];
        ax.YTick = [];
        temp_ticklengthval = [0 0];
        temp_tickdirval = 'out';
        formatHeatmapTrajPlot(ax, temp_tickdirval, temp_ticklengthval, heatmap_to_use, axesvis, rulersvis, titlesize);
    else
        ax.YTick = [];
        ax.XTick = [0 1 2 3 4];
        axesvis = true;
        rulersvis = true;
        temp_ticklengthval = [0.01 0.001];
        temp_tickdirval = 'out';
        ax.XLabel.String = 'Distance from Center of Plate (cm)';
        formatHeatmapTrajPlot(ax, temp_tickdirval, temp_ticklengthval, heatmap_to_use, axesvis, rulersvis, titlesize);
        
    end
end


% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
set(gcf, 'Position', [515   196   466   589]);
savefigname = "SuppFig_WT_Trajs_Heatmaps.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
pause(0.1);
close(gcf);

%% SUPP Figure : GFP Colony Radii

% scatter all gfps


figure(1);
clf('reset');
hold on;
count = 0;
for i = 1:length(allExptData)
    if i == 16 || i == 24
        continue;
    end
    if ~isempty(allExptData(i).scanfolder)
        for j = 1:length(allExptData(i).imfiles)
            if strcmpi(allExptData(i).genes{j}, 'gfp')
                tempiptg = allExptData(i).iptgs(j);
                if tempiptg==0
                    tempiptg= 0.01;
                end
                plot(tempiptg, allExptData(i).colrads_cm_redo(j), 'o', ...
                    'MarkerSize', 10, 'MarkerFaceColor', linecolors(2, :), ...
                    'MarkerEdgeColor', 'none');
                count = count+1;
                if allExptData(i).colrads_cm_redo(j) < 2.7
                    fprintf('Expt %g image %g\n', i, j);
                end
                    
            end
        end
    end
end
ylim([0 4.6]);
xlim([0.009 10]);
set(gca, 'XScale', 'log');
xticks([0.01 0.1 1 10]);
xticklabels({'0', '0.1', '1', '10'});
title('Colony Radii of pLacGFP vs IPTG');
xlabel('IPTG (mM)');
ylabel('Radius (cm)');

ax = gca;
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);


% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "SuppFig_All_GFP_Radii.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
pause(0.1);
close(gcf);

%% Single Gene AUCs
cd(singleExptDir);
singleGeneAUCData = loadLatestAUCData();
%% Fig 2: AUCs bar plot NOT IN USE
% 
% 
% % replace later with code to use from updated AUCs??
% % AUCs = [0.701, 0.615, 0.921, 0.879, 0.943, 0.984, 0.831];
% % gene_list = {'wt', 'gfp', 'chew', 'lrp', 'flia', 'umod', 'flgm'};
% gene_list_inds = [1 2 7 3 6 4 5];
% gene_list_order = {gene_list{gene_list_inds}};
% AUCs = zeros(1, 7);
% for i = 1:7
%     tempgene = singleGeneAUCData(i).gene;
%     gene_ind = find(strcmpi(tempgene, gene_list_order));
%     AUCs(gene_ind) = singleGeneAUCData(i).max_auc;
% end
% x = categorical(gene_list_order);
% x = reordercats(x, gene_list_order);
% 
% 
% figure(2);
% clf('reset');
% % b = bar(x, spds, 'EdgeColor', 'none');
% % title('Swarm Speed');
% 
% 
% for i = 1:7
%     linecolors_order(i, :) = linecolors(gene_list_inds(i), :);
% end
% 
% % try individual bars
% h = coloredbar2(AUCs, gene_list_order,...
%     mat2cell(linecolors_order, ones(1, 7)));
% % h = coloredbar(AUCs, gene_list_order,...
% %     mat2cell(linecolors_order, ones(1, 7)));
% ylim ([0.3, 1]);
% ylabel('AUC');
% xlabel('Strain');
% title('AUCs of Multinomial Regression Models');
% ax = gca;
% 
% % ax.FontSize = 18;
% % ax.TickLength = [0 0];
% % ax.TitleFontWeight = 'bold';
% 
% ticklength = [0 0];
% formatBarPlot(h, ax, ticklength);
% 
% % Save the figure
% set(gcf, 'Color', 'none');
% set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% savefigname = "Fig2d_AllAUCs.svg";
% export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% pause(0.1);
% close(gcf);

%% Fig 2: AUCs Dot plot instead of bar


    
% replace later with code to use from updated AUCs??
% AUCs = [0.701, 0.615, 0.921, 0.879, 0.943, 0.984, 0.831];
% gene_list = {'wt', 'gfp', 'chew', 'lrp', 'flia', 'umod', 'flgm'};
gene_list_inds = [1 2 7 3 6 5 4];
gene_list_order = {gene_list{gene_list_inds}};
AUCs = zeros(1, 7);
for i = 1:7
    tempgene = singleGeneAUCData(i).gene;
    gene_ind = find(strcmpi(tempgene, gene_list_order));
    AUCs(gene_ind) = singleGeneAUCData(i).max_auc;
end

figure(2);
clf('reset');
hold on;

for i = 1:7
    linecolors_order(i, :) = linecolors(gene_list_inds(i), :);
end

% Now plot as errorbars

for i = 1:7
    aucval = AUCs(i);
%     ind = gene_list_inds(i);
    plot(i, aucval, 'Color', linecolors_order(i, :), ...
        'DisplayName', gene_list_order{i}, ...
        'MarkerSize', 30, 'MarkerFaceColor', linecolors_order(i, :), ...
        'MarkerEdgeColor', 'none', 'Marker', 'diamond');    
%     legend_string{end+1} = gene_list_order{i};    
end
xlim([0, 8]);
xticks(1:7);
xticklabels(gene_list_order);

ylim ([0.4, 1]);
yticks([0.4 0.5 0.6 0.7 0.8 0.9 1]);
yticklabels({'0.4', '', '0.6', '', '0.8', '', '1'});
ylabel('AUC');
xlabel('Strain');
title('AUCs of Multinomial Regression Models');
ax = gca;

% ax.FontSize = 18;
% ax.TickLength = [0 0];
% ax.TitleFontWeight = 'bold';


axfontsize = 24; 
smallticks = false; 
plotgrid = false;
keeplinewidths = true;
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);


% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "Fig2d_AllAUCs.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, strrep(savefigname, '.svg', '.eps')));
pause(0.1);
close(gcf);
    
%% Fig 2: All CMs of mnr models

classLabels = categorical({'0-0.9 mM', '1-5 mM', '5-10 mM'});

for i = 1:7
    gene_curr = gene_list{i};
    tempauc = singleGeneAUCData(i).max_auc;
    figure(i);
    clf('reset');
    y_ord = singleGeneAUCData(i).y_ord;
    est_cats = singleGeneAUCData(i).est_cats;
    cm = confusionchart(double(y_ord), double(est_cats));
%     titletext = strcat(gene_curr, " Prediction Accuracy binning into ", strjoin(string(classLabels), ", "));
%     title({titletext, sprintf("mAUC: %g", tempauc), gene_curr});
    title(gene_curr);
    ax = gca;
    ax.FontSize = 14;
%     waitforbuttonpress;
    
    formatConfusMatrix(cm, cm_heatmap_to_use);

    set(gcf, 'Color', 'none');
    set(gcf, 'WindowStyle', 'normal'); % in case it was docked
    savefigname = strcat("SuppFig_AllCMs_", gene_curr, ".svg");
    export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
    pause(0.1);
    close(gcf);
    
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig 3: Load Data
cd(timelapseDir);
allExptData = loadLatestTimelapseData();

% Establish iptgs & genes to plot
iptgs = [0, 1, 5, 10];
gene_list = {'wt', 'gfp', 'chew', 'lrp', 'flia', 'umod', 'flgm'};



%% Fig 3: Combined GFP & umoD Heatmap Formatted

% umoD heatmap
load('/Users/anjali/Dropbox/Anjali_Lab_Unshared/MATLAB/umod_colormap_to_use.mat', 'umod_colormap_to_use');
figure(4);
clf('reset');
time_vec = allExptData(1).times;
%     all_dists = all_data(i).dists; 

for j = 1:5
    if j ~= 1 & j ~= 5
        continue;
    end
    ax = nexttile;

    % get averages
    avgd = allExptData(1).avgd{j};
    avgd = [flip(avgd); avgd];

    dist_vec = allExptData(1).dists{j};
    dist_vec = [flip(dist_vec) -dist_vec];

    % Cut off at 30 h
    cutoff = find(time_vec>30, 1);
    
    imagesc(time_vec(1:cutoff), dist_vec, avgd(:, 1:cutoff));
    if j == 1
        colormap(heatmap_to_use);
        caxis([0.6307    0.8457]);
    else
        colormap(umod_colormap_to_use);
        caxis([0.6757    0.8343]);
    end
    %Set the colormap once
%     colormap(colors_to_use);
%     caxis([0.63, 0.85]);
%     caxis([0.7 0.83]);
    
%     if j == 5
        cbar = colorbar;
        cbar.TickLength = 0;  cbar.AxisLocation = "out";
        drawnow;
        pause(0.1);
        cbar.Ruler.Axle.Visible = 'off';
        cbar.Box = 'off';
%     end
    
    % Make the plot pretty
    plate_name = allExptData(1).plates{j};
    plate_name = strrep(plate_name, '_', ' ');
    
    
    axesvis = true;
    rulersvis = true;
    temp_ticklengthval = [0 0];
    temp_tickdirval = 'out';
    
    if j ==5
        
        xlabel('Time elapsed (h)', 'FontSize', 15);
        ylabel('Dist from Plate Center (cm)', 'FontSize', 15);
    end
    
    formatHeatmapTrajPlot(ax, temp_tickdirval, temp_ticklengthval, heatmap_to_use, axesvis, rulersvis, titlesize);

end
% title(t, all_data(1).expt, 'FontSize', 18, 'FontWeight', 'bold');



set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% set(gcf, 'Position', [251   203   947   605]); % if on laptop
set(gcf, 'Position', [680   477   800   501]); % if on monitor
savefigname = "Fig3a_GFP_umoD_Heatmap.svg";
drawnow;
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
pause(0.1);
close(gcf);


%% Fig 3: Swarm Traj Plot


for i = 1:1 %length(allExptData)
    figure(i+5);
    clf('reset');
    hold on;
    for j = 1:length(allExptData(i).plates)
        temptraj = allExptData(i).trajs{j};
        if iscell(temptraj)
            temptraj(cellfun(@isempty, temptraj)) = {0};
            temptraj = cell2mat(temptraj);
        end
        temptraj(temptraj==0) = 1;
        if i>=9
           tempdist = allExptData(i).dists{j};
           temptraj = tempdist(round(temptraj));
        end
        temptimes = allExptData(i).times;
        tempgeneind = find(strcmpi(allExptData(i).plates{j}, gene_list));
        tempcolor = linecolors(tempgeneind, :);
        plot(temptimes, temptraj, 'Color', tempcolor, ...
            'LineWidth', 2, ...
            'DisplayName', allExptData(i).plates{j});
    end
    title(strcat(allExptData(i).expt, " ", string(allExptData(i).date)), 'Interpreter', 'none');
    xlabel('Time (h)');
    ylabel('Dist (cm)');
    axis square;
    l = legend('Location', 'northeastoutside');
    pause(1);
end

title('Swarm Trajectories');
ax = gca; 
xlim([0 40]);
% ax.YTick = [0 1 2 3 4];
% ax.XTick = [0 10 20 30 40];
keeplinewidths = false;
axfontsize = 24; 
smallticks = false; 

plotgrid = true;

formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
formatLegendText(l, gene_list, linecolors);


% ax.FontSize = 14; ax.TitleFontWeight = 'bold';


% % Save the figure
% set(gcf, 'Color', 'none');
% set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% drawnow;
% axis square;
% savefigname = "Fig3b_Trajs.svg";
% export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% pause(0.1);
% close(gcf);

%% Fig 3: Swarm Speeds (make bar plot from timelapse 1)
% x = categorical(allExptData(1).plates);
% x = reordercats(x, allExptData(1).plates);
% spds = cellfun(@mean, allExptData(1).swspds);
% % Change speedds to mm/hour--x10
% spds = spds*10;
% newlist = {'gfp', 'flia', 'chew', 'flgm', 'lrp', 'umod'};
% neworder = [1, 4, 2, 6, 3, 5];
% spds = spds(neworder);
% neworder = neworder+1;
% figure(2);
% clf('reset');
% % b = bar(x, spds, 'EdgeColor', 'none');
% % title('Swarm Speed');
% 
% 
% % h = coloredbar(spds, gene_list(2:end),...
% %     mat2cell(linecolors(2:end, :), ones(1, 6)));
% h = coloredbar2(spds, newlist,...
%     mat2cell(linecolors(neworder, :), ones(1, 6)));
% 
% ylim ([3, 9]);
% xlabel('Strain');
% title('Mean Swarm Speed at 10 mM IPTG');
% ylabel('Mean Swarm Speed (mm/h)');
% ax = gca;
% 
% % ax.FontSize = 18;
% % ax.TickLength = [0 0];
% % ax.TitleFontWeight = 'bold';
% 
% ticklength = [0 0];
% formatBarPlot(h, ax, ticklength);
% 
% % % Save the figure
% % set(gcf, 'Color', 'none');
% % set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% % savefigname = "Fig3c_Swarm_speeds.svg";
% % export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% % pause(0.1);
% % close(gcf);

%% Fig 3: Swarm Speeds DOT plot not bar
x = categorical(allExptData(1).plates);
x = reordercats(x, allExptData(1).plates);
spds = cellfun(@mean, allExptData(1).swspds);
% Change speedds to mm/hour--x10
spds = spds*10;
newlist = {'gfp', 'flia', 'chew', 'flgm', 'lrp', 'umod'};
neworder = [1, 4, 2, 6, 3, 5];
spds = spds(neworder);
neworder = neworder+1;
figure(2);
clf('reset');
hold on;

for i = 1:6
    ind = neworder(i);
    spd = spds(i); 
    plot(i, spd, 'Color', linecolors(ind, :), ...
        'DisplayName', newlist{i}, ...
        'MarkerSize', 30, 'MarkerFaceColor', linecolors(ind, :), ...
        'MarkerEdgeColor', 'none', 'Marker', 'diamond'); 
    
end


ylim ([2, 8]);
yticks([2 4 6 8]);
xlim([0, 7]);
xticks(1:6);
xticklabels(newlist);

xlabel('Strain');
% title('Mean Swarm Speed at 10 mM IPTG');
ylabel('Mean Swarm Speed (mm/h)');
ax = gca;

keeplinewidths = true;
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
set(gcf, 'Position', [698   364   709   493]);
% axis square;
drawnow;
savefigname = "Fig3c_Swarm_speeds.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
pause(0.1);
close(gcf);

%% Fig 3: Lag Cons Sw times condensed in horizontal bar plot

% Plot bar plot in certain order

spds = cellfun(@mean, allExptData(1).swspds);
swtimes = cellfun(@mean, allExptData(1).swtimes);
constimes = cellfun(@mean, allExptData(1).constimes);
lagtimes = allExptData(1).lagtime;
plate_genes = allExptData(1).plates;

% Use this to reorder it
% gene_list = {'wt', 'gfp', 'chew', 'lrp', 'flia', 'umod', 'flgm'};
% gene_list_inds = [2 6 3 7 4 5];
gene_list_inds = [5 4 7 3 6 2];
gene_list_order = {gene_list{gene_list_inds}};

y_data = zeros(6, 3);
for i = 1:6
    tempplate = plate_genes{i};
    tempind = find(strcmpi(tempplate, gene_list_order));
    disp(tempind);
    y_data(tempind, 1) = -lagtimes(i);
    y_data(tempind, 2) = constimes(i);
    y_data(tempind, 3) = swtimes(i);
end
x = categorical(gene_list_order);
x = reordercats(x, gene_list_order);
figure(3);
clf('reset');
b = barh(x, y_data,'stacked');

% change colors
for i = 1:3
    b(i).FaceColor = 'flat';
end
for i = 1:6
    genecolor = linecolors(gene_list_inds(i), :);
    b(1).CData(i, :) = genecolor/3;
    b(2).CData(i, :) = genecolor/2;
    b(3).CData(i, :) = genecolor;
end

ax = gca;
formatBarPlot(b, ax, [0 0]);
% ax.XTickLabel = {10, 5, 0, 5, 10};
ax.XLabel.String = 'Lag Time; Cons Time + Swarm Time (h)';

set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
pause(0.1);
set(gcf, 'Position', [350 214 873 466]);
% axis square;
drawnow;

savefigname = "Fig3d_Horiz_Bar.svg";
% export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% export_fig(gcf, fullfile(savefigdir, "Fig3d_Horiz_Bar.eps"));
% saveas(gcf,'/Users/anjali/Dropbox/Anjali_Updates/_Swarming_Manuscript/_High_res_figure_copies/Individual_Panels/Fig3d_Horiz_Bar.svg')
% pause(0.1);
% set(gcf, 'Color', 'white');
% set(gcf, 'WindowStyle', 'docked');

%% Fig 3: Lag Times
plotfld = 'lagtime';
phase_subset = 'all';
lagvals = getMsmtsToPlot(allExptData, plotfld, iptgs, gene_list, phase_subset);
legend_string = {};
% Decide order to plot in: gfp, umoD, fliA, cheW, flgm, lrp
order_inds = [2, 6, 5, 3, 7, 4];
figure(8);
clf('reset');
hold on;
for i = 2:7
    ind = order_inds(i-1);
%     plot(repelem(i, length(lagvals{ind, 4})), lagvals{ind, 4}, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'none', ...
%             'MarkerEdgeColor', linecolors(ind, :));
%     plot(i, mean(lagvals{ind, 4}), 'd', 'MarkerSize', 12, 'MarkerFaceColor', linecolors(ind, :), ...
%             'MarkerEdgeColor', 'none');
    meanval = mean(lagvals{ind, 4});
    stdval = std(lagvals{ind, 4});
    semval = stdval/sqrt(length(lagvals{ind, 4}));
    errorbar(i, meanval, semval, '--', 'Color', linecolors(ind, :), ...
        'DisplayName', gene_list{i}, 'LineWidth', 3, 'CapSize', 15, ...
        'MarkerSize', 30, 'MarkerFaceColor', linecolors(ind, :), ...
        'MarkerEdgeColor', 'none', 'Marker', 'diamond');    
    legend_string{end+1} = gene_list{ind};    
end
xlim([1, 8]);

ylabel(strcat('Lag Time', ' (h)'), 'FontSize', 24);
xticks(2:7);
xticklabels(gene_list(order_inds));
yticks([0 4 8 12]);
ax = gca;
ylim([0, 14]);
% 
% legend(legend_string{:});
% % ax.Legend.String = legend_string;
% ax.Legend.Location = 'northeastoutside';

keeplinewidths = true;
axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked

set(gcf, 'Position', [698   364   709   493]);
% axis square;
drawnow;

savefigname = "Fig3d_Lag_times.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
pause(0.1);
close(gcf);

%% Fig 3: Cons Times
plotfld = 'constimes';
phase_subset = 'all';
consvals = getMsmtsToPlot(allExptData, plotfld, iptgs, gene_list, phase_subset);

% Decide order to plot in: gfp, umoD, lrp, flgm, chew, flia
order_inds = [2, 7, 4, 6, 3, 5];
figure(9);
clf('reset');
hold on;
legend_string = {};
for i = 2:7
    ind = order_inds(i-1);
    disp(ind);
%     plot(repelem(i, length(consvals{ind, 4})), consvals{ind, 4}, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'none', ...
%             'MarkerEdgeColor', linecolors(ind, :));
%     legend_string{end+1} = '';
    meanval = mean(consvals{ind, 4});
    stdval = std(consvals{ind, 4});
    semval = stdval/sqrt(length(consvals{ind, 4}));
    errorbar(i, meanval, semval, '--', 'Color', linecolors(ind, :), ...
        'DisplayName', gene_list{i}, 'LineWidth', 3, 'CapSize', 15, ...
        'MarkerSize', 30, 'MarkerFaceColor', linecolors(ind, :), ...
        'MarkerEdgeColor', 'none', 'Marker', 'diamond');
%     plot(i, mean(consvals{ind, 4}), 'd', 'MarkerSize', 12, 'MarkerFaceColor', linecolors(ind, :), ...
%             'MarkerEdgeColor', 'none');
    legend_string{end+1} = gene_list{ind};
end
xlim([1, 8]);
ylabel('Consolidation Time (h)', 'FontSize', 24);
% title(strcat(plotfld, ' All Genes (10 mM IPTG)'), 'FontSize', 24);
xticks(2:7);
xticklabels(gene_list(order_inds));
ax = gca;
set(ax, 'FontSize', 14);
ylim([1, 5]);
yticks([1 2 3 4 5]);

% legend(legend_string{:});
% % ax.Legend.String = legend_string;
% ax.Legend.Location = 'northeastoutside';

axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "Fig3e_Cons_times.svg";
set(gcf, 'Position', [698   364   709   493]);
% axis square;
drawnow;

export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
pause(0.1);
close(gcf);

%% To figure out swarm times for rev response

plotfld = 'constimes';
phase_subset = 'all';
skip_expts = [8, 9];
consvals = getMsmtsToPlot(allExptData, plotfld, iptgs, gene_list, phase_subset, skip_expts);
swvals = getMsmtsToPlot(allExptData, 'swtimes', iptgs, gene_list, phase_subset, skip_expts);
sumvals = consvals;

for rownum = 1:size(sumvals, 1)
    for colnum = 1:size(sumvals, 2)
        temp_cons = consvals{rownum, colnum};
        temp_sw = swvals{rownum, colnum};
        temp_sum = temp_sw;
        for j = 1:min(length(temp_sw), length(temp_cons))
            temp_sum(j) = temp_sw(j) + temp_cons(j);
        end
        sumvals{rownum, colnum} = temp_sum;

    end
end
consvals = sumvals;
% Decide order to plot in: gfp, umoD, lrp, flgm, chew, flia
order_inds = [2, 7, 4, 6, 3, 5];
figure(9);
clf('reset');
hold on;
legend_string = {};
for i = 2:7
    ind = order_inds(i-1);
    disp(ind);
    plot(repelem(i, length(consvals{ind, 4})), consvals{ind, 4}, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', linecolors(ind, :));
    legend_string{end+1} = '';
    meanval = mean(consvals{ind, 4});
    stdval = std(consvals{ind, 4});
    semval = stdval/sqrt(length(consvals{ind, 4}));
    errorbar(i, meanval, semval, '--', 'Color', linecolors(ind, :), ...
        'DisplayName', gene_list{i}, 'LineWidth', 3, 'CapSize', 15, ...
        'MarkerSize', 30, 'MarkerFaceColor', linecolors(ind, :), ...
        'MarkerEdgeColor', 'none', 'Marker', 'diamond');
%     plot(i, mean(consvals{ind, 4}), 'd', 'MarkerSize', 12, 'MarkerFaceColor', linecolors(ind, :), ...
%             'MarkerEdgeColor', 'none');
    legend_string{end+1} = gene_list{ind};
end
xlim([1, 8]);
ylabel('Consolidation Time (h)', 'FontSize', 24);
% title(strcat(plotfld, ' All Genes (10 mM IPTG)'), 'FontSize', 24);
xticks(2:7);
xticklabels(gene_list(order_inds));
ax = gca;
set(ax, 'FontSize', 14);
% ylim([1, 5]);
yticks(1:10);
plotgrid = true;
% legend(legend_string{:});
% % ax.Legend.String = legend_string;
% ax.Legend.Location = 'northeastoutside';

axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);

%% Fig 3: Swarm Times
plotfld = 'swtimes';
phase_subset = 'all';
swvals = getMsmtsToPlot(allExptData, plotfld, iptgs, gene_list, phase_subset);

% Decide order to plot in: gfp, flgm, umod, chew, lrp, flia
order_inds = [2, 7, 6, 3, 4, 5];
figure(10);
clf('reset');
hold on;
legend_string = {};
for i = 2:7
    ind = order_inds(i-1);
    disp(ind);
%     plot(repelem(i, length(swvals{ind, 4})), swvals{ind, 4}, 'o', 'MarkerSize', 12, 'MarkerFaceColor', 'none', ...
%             'MarkerEdgeColor', linecolors(ind, :));
% %     legend_string{end+1} = '';
%     plot(i, mean(swvals{ind, 4}), 'd', 'MarkerSize', 12, 'MarkerFaceColor', linecolors(ind, :), ...
%             'MarkerEdgeColor', 'none');
        
    
    meanval = mean(swvals{ind, 4});
    stdval = std(swvals{ind, 4});
    semval = stdval/sqrt(length(swvals{ind, 4}));
    errorbar(i, meanval, semval, '--', 'Color', linecolors(ind, :), ...
        'DisplayName', gene_list{i}, 'LineWidth', 3, 'CapSize', 15, ...
        'MarkerSize', 30, 'MarkerFaceColor', linecolors(ind, :), ...
        'MarkerEdgeColor', 'none', 'Marker', 'diamond');    
    legend_string{end+1} = gene_list{ind};
end
xlim([1, 8]);
ylabel('Swarm Time (h)', 'FontSize', 24);
% title(strcat(plotfld, ' All Genes (10 mM IPTG)'), 'FontSize', 24);
xticks(2:7);
xticklabels(gene_list(order_inds));
ax = gca;
set(ax, 'FontSize', 14);
ylim([0, 7]);
yticks([0 2 4 6]);

% legend(legend_string{:});
% % ax.Legend.String = legend_string;
% ax.Legend.Location = 'northeastoutside';

axfontsize = 24; 
smallticks = false; 
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
set(gcf, 'Position', [698   364   709   493]);
% axis square;
drawnow;
savefigname = "Fig3f_Swarm_times.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
pause(0.1);
close(gcf);

%% Timelapse SUPP HEATMAPS 

for i = 2:6
    figure(i);
    gene_curr = allExptData(i).expt;
    num_plates = length(allExptData(i).avgd);
    time_vec = allExptData(i).times;
%     all_dists = all_data(i).dists; 
    t = tiledlayout(num_plates, 1);
    
    for j = 1:num_plates
        ax = nexttile;
        
        % get averages
        avgd = allExptData(i).avgd{j};
        avgd = [flip(avgd); avgd];
        
        dist_vec = allExptData(i).dists{j};
        if iscell(dist_vec)
            dist_vec = dist_vec{end};
        end
        dist_vec = [flip(dist_vec) -dist_vec];
        
%         imagesc(time_vec, dist_vec, avgd);

        % Cut off or repeat until 50 h
        % Cut off at 50 h
        cutoff = find(time_vec>50, 1);
        
        if isempty(cutoff)
            % let's add?
            gapval = 50 - time_vec(end);
            diffval = mean(diff(time_vec));
            extendval = round(gapval/diffval);
            timevec2 = zeros(1, length(time_vec) + extendval);
            timevec2(1:length(time_vec)) = time_vec;
            for l = (length(time_vec)+1): length(timevec2)
                timevec2(l) = timevec2(l-1)+diffval;
            end
            avgd2 = zeros(size(avgd, 1), length(timevec2));
            avgd2(:, 1:(size(avgd, 2))) = avgd;
%             avgd2(:, (size(avgd, 2) + 1):end) = repmat(avgd(:, end), 1, ...
%                 length(timevec2) - (size(avgd, 2)));
            avgd2(:, (size(avgd, 2) + 1):end) = NaN;
            imAlpha=ones(size(avgd2));
            imAlpha(isnan(avgd2))=0;
            imagesc(timevec2, dist_vec, avgd2,'AlphaData',imAlpha);
        else
            imagesc(time_vec(1:cutoff), dist_vec, avgd(:, 1:cutoff));
        end


        %Set the colormap once
        colormap(heatmap_to_use);
        
        
        if j == 4
            caxis([0.5, 0.87]);
            cbar = colorbar;
            cbar.TickLength = 0; cbar.Box = 'off'; cbar.AxisLocation = "out";
            drawnow;
            pause(0.1);
            cbar.Ruler.Axle.Visible = 'off';
        end
        
        % Make the plot pretty
        plate_name = allExptData(i).plates{j};
        plate_name = strrep(plate_name, '_', ' ');
%         title(plate_name,'FontSize',14, 'FontWeight', 'normal', 'Interpreter', 'none');
%         ax.LineWidth = 1.5;
%         ax.Box = 'off';
%         xlim([0 60]);
        if j ~= 4
            ax.XTick = [];
        else
            ax.XTick = [0 10 20 30 40 50];
        end
        axesvis = true;
        rulersvis = false;
        temp_ticklengthval = [0 0];
        temp_tickdirval = 'out';
        formatHeatmapTrajPlot(ax, temp_tickdirval, temp_ticklengthval, heatmap_to_use, axesvis, rulersvis, titlesize);

    end
    title(t, allExptData(i).expt, 'FontSize', 18, 'FontWeight', 'bold');
    
    xlabel(t, 'Time elapsed (h)', 'FontSize', 15);
    ylabel(t, 'Dist from Plate Center (cm)', 'FontSize', 15);
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    
    % Save the figure
    set(gcf, 'Color', 'none');
    set(gcf, 'WindowStyle', 'normal'); % in case it was docked
    set(gcf, 'Position', [ 439   280   565   520]);
    drawnow;
    savefigname = strcat("SuppFig_TimelapseHeatmaps_", gene_curr, ".svg");
    export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
    pause(0.1);
    close(gcf);
end


%% Timelapse supp plots vs iptg

plotflds = {'lagtime', 'swtimes', 'constimes', 'swspds'};
phase_subset_options = {'all', 'first', 'middle', 'last'};
tempxvals = [0.1 1 5 10];
iptgs = [0 1 5 10];
for subset_ind = 1:4
    phase_subsets = {'all', phase_subset_options{subset_ind}, ...
        phase_subset_options{subset_ind}, phase_subset_options{subset_ind}};
% phase_subsets = {'all', 'first', 'first', 'first'};
    for plotfld_ind = 1:length(plotflds)
        plotfld = plotflds{plotfld_ind};
        phase_subset = phase_subsets{plotfld_ind};
        skip_expts = [8, 9];
        plotvals = getMsmtsToPlot(allExptData, plotfld, iptgs, gene_list, phase_subset, skip_expts);

        figure(plotfld_ind);
        clf('reset');
        hold on;
        for i = 3:7
    %         figure(i);
    %         clf('reset');
    %         hold on;
            meanvals = zeros(1, 4);
            for j = 1:length(iptgs)
                tempvals = plotvals{i, j};
                plot(repelem(tempxvals(j), length(tempvals)), tempvals, 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', linecolors(i, :), 'DisplayName', '');
                meanvals(j) = nanmean(tempvals);
            end

%             % Plot mean line through them
%             plot((1:4), meanvals, 'o-','LineWidth', 5, 'Color', linecolors(i, :), 'MarkerSize', 15, ...
%                 'MarkerFaceColor', linecolors(i, :));
            meanlines(i) = plot(tempxvals, meanvals, 'o-','LineWidth', 5, 'Color', linecolors(i, :), 'MarkerSize', 15, ...
                'MarkerFaceColor', linecolors(i, :), 'DisplayName', gene_list{i});
            title(strcat(plotfld, " ", phase_subset));
    %         title(strcat(gene_list{i}, " ", plotfld, " (h)"), 'FontSize', 24);
%             xlim([0, 4.5]);
        %     ylim([0, 10]);
            
            ax = gca;
            axfontsize = 24; 
            smallticks = false; 
            formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
                keeplinewidths, axfontsize, smallticks, titlesize);

            % take a look at the figure
            drawnow;
%             waitforbuttonpress;

    %         % Save the figure
    %         set(gcf, 'Color', 'none');
    %         set(gcf, 'WindowStyle', 'normal'); % in case it was docked
    %         set(gcf, 'Position', [698   364   709   493]);
    %         savefigname = strcat("SuppFig_", plotfld, "vsIPTG_", gene_list{i}, ".svg");
    % %             savefigname = "SuppFig_LagTimes_vsIPTG_lrp.svg";
    %         drawnow;
    %         export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
    %         pause(0.1);
    %         close(gcf);

        end
        l = legend(meanlines(3:7), 'Location', 'eastoutside');
        formatLegendText(l, gene_list, linecolors);
        ax.XLim = [0.09 10];
        ax.XTick = [0.1 1 5 10];
        ax.XTickLabel = {'0', '1', '5', '10'};
        
        if strcmpi(plotfld, 'constimes')
            ax.YLim = [1 8];
        end

        set(ax, 'XScale', 'log');
        
        disp(strcat("SuppFig_", plotfld, "_vsIPTG_", phase_subset, ".svg"));
%         waitforbuttonpress;
        
         % Save the figure
        set(gcf, 'Color', 'none');
        set(gcf, 'WindowStyle', 'normal'); % in case it was docked
        set(gcf, 'Position', [698   364   709   493]);
%         savefigname = strcat("SuppFig_", plotfld, "vsIPTG_", gene_list{i}, ".svg");
%             savefigname = "SuppFig_LagTimes_vsIPTG_lrp.svg";
        drawnow;
        
        savefigname = strcat("SuppFig_", plotfld, "_vsIPTG_", phase_subset, ".svg");
        export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
        
        set(gcf, 'Color', 'white');
        drawnow;
        pause(0.5);
        savefigname = strcat("SuppFig_", plotfld, "_vsIPTG_", phase_subset, ".png");
        export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
        
        pause(0.1);
        close(gcf);
    end
    
end

%% Fig 3e LRP Swarming Phase--Front CV vs IPTG

% Plot front CV vs iptg
tempxvals = [0.1 1 5 10];
iptgs = [0 1 5 10];
phase_val = 2; %use 1 for cons
sw_dens_expts = [3 4 5 6 7];
figure(1);
clf('reset');
hold on;
for i = 7:7 %2:length(allExptData)
    if ~strcmpi(allExptData(i).expt, 'lrp')
        continue;
    end
    dens_vals = zeros(1, 4);
    dens_stds = dens_vals;
    dens_sems = dens_vals;
    % For each iptg, get the mean & the std
    for j = 1:4
        temptrajlabels = allExptData(i).traj_labels{j};
        tempdensvals = allExptData(i).front_cv_trajs{j};
        
        tempdensvals = tempdensvals(temptrajlabels==phase_val);
        
%         % Let's find the values of the first swarm phase:
%         first_cons_bound = find(temptrajlabels==phase_val, 1);
%         for k = first_cons_bound:length(temptrajlabels)
%             if temptrajlabels(k) ~= phase_val
%                 first_cons_end = k;
%                 break;
%             end
%         end
%         tempdensvals = tempdensvals(first_cons_bound:first_cons_end);
%         
        dens_vals(j) = nanmean(tempdensvals);
%         dens_vals(j) = mean(tempdensvals);
        dens_stds(j) = nanstd(tempdensvals);
        dens_sems(j) = nanstd(tempdensvals)/sqrt(length(tempdensvals));
        
    end
    tempgene = allExptData(i).expt;
    gene_ind = find(strcmpi(tempgene, gene_list));
    % Plot
%     meanlines(i) = plot(tempxvals, dens_vals, 'LineWidth', 2, 'Color', linecolors(gene_ind, :), ...
%         'DisplayName', tempgene);
    e = errorbar(tempxvals, dens_vals, dens_sems, ...
        'd-', 'MarkerSize', 30, ...
        'Color', linecolors(gene_ind, :), ...
        'MarkerEdgeColor', 'none',...
        'MarkerFaceColor', linecolors(gene_ind, :),...
        'LineWidth', 3, 'DisplayName', '');
    
    % %       'Marker', '_', 'MarkerSize', 10, 'CapSize', 20, 
    
end
% legend('Location', 'bestoutside');;
drawnow;
e.CapSize = 20;
xlim([-1 10]);
ax = gca;
ax.XTick = [0.1 1 5 10];
ax.XTickLabel = {'0', '1', '5', '10'};
xlabel('IPTG (mM)');
ylabel('Colony Front CV during Swarming Phase');
title('lrp Swarm Phase Front CV vs IPTG');
axfontsize = 24; 
smallticks = false; 
keeplinewidths = true;
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
 ax.YAxis.Exponent = 0;
% set(gca, 'FontSize', 18);
set(gcf, 'Color', 'white');
% set(gca, 'XScale', 'log');

% % For saving:
% set(gcf, 'Color', 'none');
% set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% pause(0.1);
% set(gcf, 'Position', [680   258   605   540]);
% drawnow;
% savefigname = strcat("Fig3f_SwarmFrontCVlrp.svg");
% export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% set(gcf, 'Color', 'white');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Combo Gene Data
cd(comboExptDir);
allExptData = loadLatestComboSensorData();
%% Fig 4: Traj heatmaps

iptgs = [0, 0.5, 1];
aras = [0, 0.05, 0.1];

% This plot should create the heatmap trajectories of the plates pictured in 4b. 
% Therefore, need the dates/image names of those.
% Currently in use:
% cheW_umoD_0_IPTG_0_ara_002.tif, 2-20-20
% cheW_umoD_0_IPTG_0.1_ara_002.tif, 2-20-20
% from date 1-14-21:
% 0 IPTG 0.05 ara: 002
% 0.5 IPTG 0 ara: 001
% 0.5/0.05 002
% 0.5/0.1: v2002
% 1/0: 002
% 1/0.05: 002
% 1/0.1: 002

% Get the trajs specified in the images, get the *distances* for each traj,
% then interpolate to get each traj on the range 0 to 4 at some
% predetermined interval
all_trajs = {};
all_labels = {};

% Choose expts/image indexes in expts to plot; first column is expt inds,
% 2nd column is image inds; rows are IPTG 0, 0.5, 1; columns are ara 0,
% 0.05, 0.1 %
inds = ones(9, 2); % use 1st expt for most (1-14-21)
inds([1, 3], 1) = 7; %use 2-20-20 expt for two controls

% Get the image indexes
tempfiles = allExptData(7).imfiles; %2-20-20
inds(1, 2) = find(contains(tempfiles, 'cheW_umoD_0_IPTG_0_ara_002.tif'));
inds(3, 2) = find(contains(tempfiles, 'cheW_umoD_0_IPTG_0.1_ara_002.tif'));
% filling in the rest manually
inds(2, 2) = 15;
inds(4:9, 2) = [8, 2, 5, 27, 19, 23]';
all_labels_txt = {'0i/0a', '0i/0.05a', '0I/0.1a', '0.5i/0a', '0.5i/0.05a', ...
    '0.5i/0.1a', '1i/0a', '1i/0.05a', '1i/0.1a'};
xvals = 0:0.0075:4;

% Iterate over each image
for i = 1:9
    exptind = inds(i, 1);
    imind = inds(i, 2);
    % Get the trajectory
    temptraj = allExptData(exptind).avgd{imind};
    if isnan(temptraj(1))
        temptraj(1) = temptraj(2);
    end
    temprad = allExptData(exptind).colrads_pix(imind);
    tempradcm = allExptData(exptind).colrads_cm(imind);
    % For now, fill in all values outside of colony with 1's
    if length(temptraj)==1000
        temptraj(round(temprad):end) = 1;
    else
        newtraj = ones(1, 1000);
        newtraj(1:length(temptraj)) = temptraj;
        temptraj = newtraj;
    end
    
    % Get the prior dist vector
    tempdist = allExptData(exptind).dists{imind};
    % Make the new trajectory via interpolation (from 0 to 4 cm)
    tempvec = interp1(tempdist, temptraj, xvals);
    % Now reset anything outside of the colony to NaN
    radind = find((xvals>tempradcm), 1);
    tempvec(radind:end) = NaN;
    
    % Unsure if needed--if the first value is NaN after interpolation need
    % this?
%     tempvec = tempvec(2:end); % for some reason interp1 makes first value NaN
%     xvals = xvals(2:end);

    % Store the new trajectory in all_trajs
    all_trajs{i} = tempvec;
end

clear exptind imind i tempvec tempradcm tempdist temptraj temprad tempfiles
all_trajs = reshape(all_trajs, [3 3]);
all_trajs = all_trajs';

% Make Labels for heatmap
iptgs = [0, 0.5, 1];
aras = [0, 0.05, 0.1];

% SET The Colormap for the whole figure

% Make Heatmaps
figure(7);
% set(gcf, 'Visible', 'off');
clf('reset');
t = tiledlayout(3, 1);

% Want to make 3 different groups (one for each IPTG concentration)

for i = 1:3
    tempara = aras(i);
    ax = nexttile;
    
    temptrajs = (vertcat(all_trajs{:, i}));
    yvals = iptgs;
    
    
    % Imagesc approach

    imAlpha=ones(size(temptrajs));
    imAlpha(isnan(temptrajs))=0;
    imagesc(xvals, yvals, temptrajs,'AlphaData',imAlpha);
    colormap(heatmap_to_use);
    caxis([0.65 1.0]);
    
    title(sprintf('%g%% Arabinose', tempara));
    
    if i == 2
        ax.YLabel.String = 'IPTG (mM)';
        c = colorbar;
        c.TickLength = 0; c.Box = 'off';
        c.AxisLocation = "out";
        c.Ticks = [0.7 0.85 1];
        c.TickLabels = {'0.7', '0.85', '1'};
        drawnow;
        pause(0.1);
        c.Ruler.Axle.Visible = 'off';
    end
    
    axesvis = true;
    rulersvis = true;
    xticks([0 1 2 3 4]);
    temp_ticklengthval = [0 0];
    temp_tickdirval = 'out';
    formatHeatmapTrajPlot(ax, temp_tickdirval, temp_ticklengthval, heatmap_to_use, axesvis, rulersvis, titlesize);
    % turn off only the y rulers
    ax.YRuler.Axle.Visible = 'off';
end
ax.XLabel.String = 'Distance From Center (cm)';

set(gcf, 'Visible', 'on');
% 
% % Save the figure
% set(gcf, 'Color', 'none');
% set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% % export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% exportgraphics(gcf, fullfile(savefigdir, 'Fig4_c_Traj_Heatmaps.eps' ));
% pause(0.1);
% close(gcf);

%% Fig 4: Col area

% Use this section to make a plot of colony plate area for desired conditions. 
% For the top 3 conditions which reach the edge, set them to some arbitrary number (maybe max radius of plate? ~4 cm)

clear all_areas

% Combos to use decided above:
% iptgs = [0, 0.5, 1];
% aras = [0, 0.05, 0.1];

iptgs = [0, 2.5, 5.0];
aras = [0, 0.1, 0.2];

% Generate matrix of combos to use
temp1 = repelem(iptgs, length(aras));
temp2 = repmat(aras', length(iptgs), 1);
use_combos = [temp1', temp2]; % Creates vertical 9x2 vector where column 1 is iptgs & column 2 is aras
all_combos = cell(9, 1);
for i = 1:length(use_combos)
    all_combos{i} = use_combos(i, :);
end
all_areas = cell(3, 3);

% Make a cell array to store rads of each plate matching condition for each
% combo; each row will be iptg, columns will be ara

% Iterate over each combo of IPTG & ara to get inds of plates for each 
for i = 1:numel(all_areas)
    iptg = all_combos{i}(1);
    ara = all_combos{i}(2);
    fprintf('Combo %g iptg, %g ara; %g of %g\n', iptg, ara, i, numel(all_areas));
    if isempty(all_areas{i})
        for j = 1:10 %length(allExptData)
            % Iterate over allExptData looking for plates at that comb of conditions
            fprintf('Expt %g of %g\n', j, length(allExptData));
            % IMPORTANT: For the 1/14/21 expt, the 0 IPTG/0 ara and
            % 0/IPTG/0.1 ara plates were outliers. don't use
            if iptg==0 && (ara == 0 | ara == 0.1) && j == 1
                continue;
            end
            
            % Find if any IPTG/aras are in that combo
            tempcombos = [allExptData(j).iptgs', allExptData(j).aras'];
            inds = ((allExptData(j).iptgs==iptg) & (allExptData(j).aras==ara));
            if j == 1
                outlier_inds = contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.1_a_v2004') |...
                    contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.05_a_004');
                inds = inds & ~outlier_inds;
            end
            inds = find(inds);
            
            % If so store the colony percent areas in each
            tempareas = allExptData(j).percent_area(inds);
            tempareas(tempareas>1) = 1;
            all_areas{i} = [all_areas{i} tempareas];
            
        end % Move to next expt 
    end % end if all_rads is empty section
    
end % Move to next iptg/ara combo
all_areas = all_areas';


meanareas = cellfun(@nanmean, all_areas);
% Now, generate heatmap from the mean rad cms
figure(1);

clf('reset');

h = heatmap(aras, iptgs, meanareas, 'GridVisible', 'off');
% colormap(colors_to_use);
title('Mean Colony Area (% of full plate)');
ax = gca;
set(ax,'FontSize',12);
ax.XLabel = 'Arabinose (%)';
ax.YLabel = 'IPTG (mM)';
% ax.LineWidth = 1.5;
% ax.Box = 'off';

formatHeatMapChart(h, h_colormap_to_use)
set(gcf, 'Visible', 'on');

% % Save the figure
% set(gcf, 'Color', 'none');
% set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% % export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% exportgraphics(h, fullfile(savefigdir, 'Fig4d_ColRads_Heatmap.eps' ));
% pause(0.1);
% close(gcf);

%% Fig 4: Get precise n of images

% Use this section to make a plot of colony plate area for desired conditions. 
% For the top 3 conditions which reach the edge, set them to some arbitrary number (maybe max radius of plate? ~4 cm)

clear all_areas

% Combos to use decided above:
% iptgs = [0, 0.5, 1];
% aras = [0, 0.05, 0.1];

iptgs = [0, 2.5, 5.0];
aras = [0, 0.1, 0.2];

% Generate matrix of combos to use
temp1 = repelem(iptgs, length(aras));
temp2 = repmat(aras', length(iptgs), 1);
use_combos = [temp1', temp2]; % Creates vertical 9x2 vector where column 1 is iptgs & column 2 is aras
all_combos = cell(9, 1);
for i = 1:length(use_combos)
    all_combos{i} = use_combos(i, :);
end
all_areas = cell(3, 3);

% Make a cell array to store rads of each plate matching condition for each
% combo; each row will be iptg, columns will be ara

% Iterate over each combo of IPTG & ara to get inds of plates for each 
for i = 1:numel(all_areas)
    iptg = all_combos{i}(1);
    ara = all_combos{i}(2);
    fprintf('Combo %g iptg, %g ara; %g of %g\n', iptg, ara, i, numel(all_areas));
    count = 0;
    if isempty(all_areas{i})
        for j = 1:10 %length(allExptData)
            % Iterate over allExptData looking for plates at that comb of conditions
%             fprintf('Expt %g of %g\n', j, length(allExptData));
            % IMPORTANT: For the 1/14/21 expt, the 0 IPTG/0 ara and
            % 0/IPTG/0.1 ara plates were outliers. don't use
            if iptg==0 && (ara == 0 | ara == 0.1) && j == 1
                continue;
            end
            
            % Find if any IPTG/aras are in that combo
            tempcombos = [allExptData(j).iptgs', allExptData(j).aras'];
            inds = ((allExptData(j).iptgs==iptg) & (allExptData(j).aras==ara));
            if j == 1
                outlier_inds = contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.1_a_v2004') |...
                    contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.05_a_004');
                inds = inds & ~outlier_inds;
            end
            inds = find(inds);
            
            % If so store the colony percent areas in each
            tempareas = allExptData(j).percent_area(inds);
            tempareas(tempareas>1) = 1;
            all_areas{i} = [all_areas{i} tempareas];
            count = count + length(tempareas);
        end % Move to next expt 
    end % end if all_rads is empty section
    disp(count);
end % Move to next iptg/ara combo
all_areas = all_areas';
% Note the output:
% Combo 0 iptg, 0 ara; 1 of 9: 41
% Combo 0 iptg, 0.1 ara; 2 of 9: 23
% Combo 0 iptg, 0.2 ara; 3 of 9: 27
% Combo 2.5 iptg, 0 ara; 4 of 9: 35
% Combo 2.5 iptg, 0.1 ara; 5 of 9: 22
% Combo 2.5 iptg, 0.2 ara; 6 of 9: 14
% Combo 5 iptg, 0 ara; 7 of 9: 31
% Combo 5 iptg, 0.1 ara; 8 of 9: 23
% Combo 5 iptg, 0.2 ara; 9 of 9: 26


%% Fig 4 SUPP Col Area with SEMs
% Let's make a bar plot with groups of 3 bars, the SEMs, and the individual
% data scattered over it?

% get the values
clear all_areas

% Combos to use decided above:
% iptgs = [0, 0.5, 1];
% aras = [0, 0.05, 0.1];

iptgs = [0, 2.5, 5.0];
aras = [0, 0.1, 0.2];

% Generate matrix of combos to use
temp1 = repelem(iptgs, length(aras));
temp2 = repmat(aras', length(iptgs), 1);
use_combos = [temp1', temp2]; % Creates vertical 9x2 vector where column 1 is iptgs & column 2 is aras
all_combos = cell(9, 1);
for i = 1:length(use_combos)
    all_combos{i} = use_combos(i, :);
end
all_areas = cell(3, 3);

% Make a cell array to store rads of each plate matching condition for each
% combo; each row will be iptg, columns will be ara

% Iterate over each combo of IPTG & ara to get inds of plates for each 
for i = 1:numel(all_areas)
    iptg = all_combos{i}(1);
    ara = all_combos{i}(2);
%     fprintf('Combo %g iptg, %g ara; %g of %g\n', iptg, ara, i, numel(all_areas));
    if isempty(all_areas{i})
        for j = 1:10 %length(allExptData)
            % Iterate over allExptData looking for plates at that comb of conditions
            
%             fprintf('Expt %g of %g\n', j, length(allExptData));
            % IMPORTANT: For the 1/14/21 expt, the 0 IPTG/0 ara and
            % 0/IPTG/0.1 ara plates were outliers. don't use
            if iptg==0 && (ara == 0 | ara == 0.1) && j == 1
                continue;
            end
            
            % Find if any IPTG/aras are in that combo
            tempcombos = [allExptData(j).iptgs', allExptData(j).aras'];
            inds = ((allExptData(j).iptgs==iptg) & (allExptData(j).aras==ara));
            if j == 1
                outlier_inds = contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.1_a_v2004') |...
                    contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.05_a_004');
                inds = inds & ~outlier_inds;
            end
            inds = find(inds);
            
            % If so store the colony percent areas in each
            tempareas = allExptData(j).percent_area(inds);
            tempareas(tempareas>1) = 1;
            all_areas{i} = [all_areas{i} tempareas];
            
        end % Move to next expt 
    end % end if all_rads is empty section
    
end % Move to next iptg/ara combo
all_areas = all_areas';

% Get the means and sems
meanareas = cellfun(@nanmean, all_areas);
stdareas = cellfun(@std, all_areas);
semareas = stdareas;
for i = 1:numel(all_areas)
    semareas(i) = stdareas(i)./sqrt(length(all_areas{i}));
end

% Bar plot?
figure(2);
clf('reset');
b = bar(meanareas);
% Format
% Add the scatter of sems over it

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(meanareas);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% % Plot the errorbars
% e = errorbar(x',meanareas,semareas,'k','linestyle','none');
% 
% 
% % FORMATT whoo
% for i = 1:length(e)
%     e(i).LineWidth = 2;
%     e(i).CapSize = 5;
% end
ax = gca;
ticklength = [0 0];
ax.XTickLabel = {'0 IPTG', '2.5 IPTG', '5 IPTG'};
title('Mean Areas w SEMs');

formatBarPlot(b, ax, ticklength);
% tempcolors = [8, 48, 107; 85, 125, 184; 161, 201, 255]/255;
tempcolors = [53, 75, 124; 85, 125, 184; 161, 201, 255]/255;

for i = 1:3
    b(i).FaceColor = 'flat'; 
end

for i = 1:3
    b(i).CData(1, :) = tempcolors(i, :);
    b(i).CData(2, :) = tempcolors(i, :);
    b(i).CData(3, :) = tempcolors(i, :);
end
% for i = 1:3
%     b(i).FaceColor = 'none'; 
% end
% for i = 1:3
%     
%     b(i).EdgeColor = 'flat';
%     b(i).LineWidth = 2; 
% end

% b.DisplayName


% Fig 4 SUPP Col Area with SEMs: Add SCATTER 01-24-23


% Scatter all the data points
x_transpose = x';
for i = 1:numel(all_areas)
    tempx = x_transpose(i);
    temp_areas = all_areas{i};
    s(i) = scatter(repelem(tempx, length(temp_areas)), ...
            temp_areas, 25, 'filled', 'MarkerFaceColor', 'k', ... %[0.01 0.01 0.01], ...
            'MarkerFaceAlpha', 0.7, 'XJitter', 'rand', 'XJitterWidth', 0.1);
end
% Add the error bars
% Plot the errorbars
e = errorbar(x',meanareas,semareas,'k','linestyle','none');

% FORMATT whoo
for i = 1:length(e)
    e(i).LineWidth = 2;
    e(i).CapSize = 5;
end

ax = gca;
ticklength = [0 0];
ax.XTickLabel = {'0 IPTG', '2.5 IPTG', '5 IPTG'};
title('Colony Areas');
ax.YLim = [0 1.1];
% formatBarPlot(b, ax, ticklength);
% b.DisplayName

%% Try Violin plot of Areas instead

% Run the previous section first
figure(3);
clf('reset');
% % hold on;
% Scatter all the data points
for i = 1:numel(all_areas)
    tempx = i;
    temp_areas = all_areas{i};
    boxoff = true;
    al_goodplot_ad(temp_areas', i, 0.5,[0.5 0.5 1],'bilateral',[],std(temp_areas(:))/1000, 0.05, boxoff);
%     s(i) = scatter(repelem(tempx, length(temp_areas)), ...
%             temp_areas, 50, 'filled', 'MarkerFaceColor', [0.3 0.5 0.8], ...
%             'MarkerFaceAlpha', 0.7, 'XJitter', 'rand', 'XJitterWidth', 0.1);
end

ylim([0 1.1]);
ax = gca;
ticklength = [0 0];
% ax.XTickLabel = {'0 IPTG', '2.5 IPTG', '5 IPTG'};
title('Mean Areas w SEMs');
% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, 'FigS9b_areas_violin.eps' ));
pause(0.1);
set(gcf, 'WindowStyle', 'docked'); % in case it was docked
set(gcf, 'Color', 'white');
disp('done');

%% Fig 4: CV heatmaps
iptgs = [0, 2.5, 5.0];
aras = [0, 0.1, 0.2];

% Generate matrix of combos to use
temp1 = repelem(iptgs, length(aras));
temp2 = repmat(aras', length(iptgs), 1);
use_combos = [temp1', temp2]; % Creates vertical 9x2 vector where column 1 is iptgs & column 2 is aras
all_combos = cell(9, 1);
for i = 1:length(use_combos)
    all_combos{i} = use_combos(i, :);
end
all_CVs = cell(3, 3);

% Make a cell array to store rads of each plate matching condition for each
% combo; each row will be iptg, columns will be ara

% Iterate over each combo of IPTG & ara to get inds of plates for each 
for i = 1:numel(all_CVs)
    iptg = all_combos{i}(1);
    ara = all_combos{i}(2);
    fprintf('Combo %g iptg, %g ara; %g of %g\n', iptg, ara, i, numel(all_CVs));
    if isempty(all_CVs{i})
        for j = 1:10 %length(allExptData)
            % Iterate over allExptData looking for plates at that comb of conditions
            fprintf('Expt %g of %g\n', j, length(allExptData));
            % IMPORTANT: For the 1/14/21 expt, the 0 IPTG/0 ara and
            % 0/IPTG/0.1 ara plates were outliers. don't use
            if iptg==0 && (ara == 0 | ara == 0.1) && j == 1
                continue;
            end
            
            % Find if any IPTG/aras are in that combo
            tempcombos = [allExptData(j).iptgs', allExptData(j).aras'];
            inds = ((allExptData(j).iptgs==iptg) & (allExptData(j).aras==ara));
            if j == 1
                outlier_inds = contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.1_a_v2004') |...
                    contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.05_a_004');
                inds = inds & ~outlier_inds;
            end
            inds = find(inds);
            
            % If so store the colony percent areas in each
            tempCVs = allExptData(j).mean_cv(inds);
            
            all_CVs{i} = [all_CVs{i} tempCVs];
            
        end % Move to next expt 
    end % end if all_rads is empty section
    
end % Move to next iptg/ara combo
all_CVs = all_CVs';


meanCVs = cellfun(@nanmean, all_CVs);

% Now, generate heatmap from the mean rad cms
figure(1);
clf('reset');
h = heatmap(aras, iptgs, meanCVs, 'GridVisible', 'off');
% colormap(colors_to_use);
title('Mean CV (Sliding Window)');
ax = gca;
set(ax,'FontSize',12);
ax.XLabel = 'Arabinose (%)';
ax.YLabel = 'IPTG (mM)';
formatHeatMapChart(h, h_colormap_to_use);

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(h, fullfile(savefigdir, 'Fig4e_CVs_Heatmap.eps' ));
pause(0.1);
close(gcf);

%% Fig 4 SUPP Figure Bar plot of CVs with SEMs

% Let's make a bar plot with groups of 3 bars, the SEMs, and the individual
% data scattered over it?

% get the values
clear all_cvs

% Combos to use decided above:
% iptgs = [0, 0.5, 1];
% aras = [0, 0.05, 0.1];

iptgs = [0, 2.5, 5.0];
aras = [0, 0.1, 0.2];

% Generate matrix of combos to use
temp1 = repelem(iptgs, length(aras));
temp2 = repmat(aras', length(iptgs), 1);
use_combos = [temp1', temp2]; % Creates vertical 9x2 vector where column 1 is iptgs & column 2 is aras
all_combos = cell(9, 1);
for i = 1:length(use_combos)
    all_combos{i} = use_combos(i, :);
end
all_cvs = cell(3, 3);

% Make a cell array to store rads of each plate matching condition for each
% combo; each row will be iptg, columns will be ara

% Iterate over each combo of IPTG & ara to get inds of plates for each 
for i = 1:numel(all_cvs)
    iptg = all_combos{i}(1);
    ara = all_combos{i}(2);
%     fprintf('Combo %g iptg, %g ara; %g of %g\n', iptg, ara, i, numel(all_areas));
    if isempty(all_cvs{i})
        for j = 1:10 %length(allExptData)
            % Iterate over allExptData looking for plates at that comb of conditions
            
%             fprintf('Expt %g of %g\n', j, length(allExptData));
            % IMPORTANT: For the 1/14/21 expt, the 0 IPTG/0 ara and
            % 0/IPTG/0.1 ara plates were outliers. don't use
            if iptg==0 && (ara == 0 | ara == 0.1) && j == 1
                continue;
            end
            
            % Find if any IPTG/aras are in that combo
            tempcombos = [allExptData(j).iptgs', allExptData(j).aras'];
            inds = ((allExptData(j).iptgs==iptg) & (allExptData(j).aras==ara));
            if j == 1
                outlier_inds = contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.1_a_v2004') |...
                    contains(allExptData(j).imfiles, 'cheWumoD_0.5_i_0.05_a_004');
                inds = inds & ~outlier_inds;
            end
            inds = find(inds);
            
            % If so store the colony percent areas in each
            tempcvs = allExptData(j).mean_cv(inds);
            if any(tempcvs==0)
                disp('found 0');
                disp(j);
                disp(inds);
            end
            all_cvs{i} = [all_cvs{i} tempcvs];
            
        end % Move to next expt 
    end % end if all_rads is empty section
    
end % Move to next iptg/ara combo
all_cvs = all_cvs';

% Get the means and sems
meancvs = cellfun(@nanmean, all_cvs);
stdcvs = cellfun(@std, all_cvs);
semcvs = stdareas;
for i = 1:numel(all_cvs)
    semcvs(i) = stdcvs(i)./sqrt(length(all_cvs{i}));
end

% Bar plot?
figure(2);
clf('reset');
b = bar(meancvs);
% Format
% Add the scatter of sems over it

hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(meancvs);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
e = errorbar(x',meancvs,semcvs,'k','linestyle','none');


% FORMATT whoo
for i = 1:length(e)
    e(i).LineWidth = 2;
    e(i).CapSize = 5;
end
ax = gca;
ax.YTick = [0, 0.01, 0.02, 0.03];
ax.XTickLabel = {'0 IPTG', '2.5 IPTG', '5 IPTG'};
ticklength = [0 0];
title('Mean CVs w SEMs');


for i = 1:3
    b(i).FaceColor = 'flat';
end
% tempcolors = [8, 48, 107; 85, 125, 184; 161, 201, 255]/255;
tempcolors = [53, 75, 124; 85, 125, 184; 161, 201, 255]/255;
for i = 1:3
    
    b(i).CData(1, :) = tempcolors(i, :);
    b(i).CData(2, :) = tempcolors(i, :);
    b(i).CData(3, :) = tempcolors(i, :);
end

formatBarPlot(b, ax, ticklength);

% for i = 1:3
%     b(i).FaceColor = 'none'; 
% end
% for i = 1:3
%     
%     b(i).EdgeColor = 'flat';
%     b(i).LineWidth = 2; 
% end
% b.DisplayName

% Bar plot of CVs with Scatter Overlaid

hold on;
% Scatter all the data points
x_transpose = x';
for i = 1:numel(all_cvs)
    tempx = x_transpose(i);
    temp_cvs = all_cvs{i};
    s(i) = scatter(repelem(tempx, length(temp_cvs)), ...
            temp_cvs, 25, 'filled', 'MarkerFaceColor', [0.01 0.01 0.01], ...
            'MarkerFaceAlpha', 0.6, 'XJitter', 'rand', 'XJitterWidth', 0.1);
end

ax.YTick = [0, 0.01, 0.02, 0.03, 0.04];
ax.YLim = [0, 0.045];
%%
% Format the colors of the barchart
for i = 1:3
    b(i).FaceColor = 'flat';
end
tempcolors = [8, 48, 107; 85, 125, 184; 161, 201, 255]/255;
for i = 1:3
    
    b(i).CData(1, :) = tempcolors(i, :);
    b(i).CData(2, :) = tempcolors(i, :);
    b(i).CData(3, :) = tempcolors(i, :);
end


%% Violin Plot of the CVs

% Run the previous section first
% Note that after the previous section all cvs is in a different order than
% all areas? fix
all_cvs = all_cvs';
figure(3);
clf('reset');
% % hold on;
% Scatter all the data points
for i = 1:numel(all_cvs)
    tempx = i;
    temp_cvs = all_cvs{i};
    boxoff = true;
    al_goodplot_ad(temp_cvs', i, 0.5,[0.5 0.5 1],'bilateral',[],std(temp_cvs(:))/1000, 0.05, boxoff);
%     s(i) = scatter(repelem(tempx, length(temp_areas)), ...
%             temp_areas, 50, 'filled', 'MarkerFaceColor', [0.3 0.5 0.8], ...
%             'MarkerFaceAlpha', 0.7, 'XJitter', 'rand', 'XJitterWidth', 0.1);
end

% ylim([0 1.1]);
ax = gca;
ticklength = [0 0];
% ax.XTickLabel = {'0 IPTG', '2.5 IPTG', '5 IPTG'};
title('Mean CVs w SEMs');
% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
% export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
exportgraphics(gcf, fullfile(savefigdir, 'FigS9b_cvs_violin.eps' ));
pause(0.1);
set(gcf, 'WindowStyle', 'docked'); % in case it was docked
set(gcf, 'Color', 'white');
disp('done');
%% Fig 4: AUCs
cd(comboExptDir);
comboGeneAUCData = loadLatestComboAUCData();
classLabels = categorical({'0-1 mM IPTG, 0-0.09% ara',...
    '1-4.9 mM IPTG, 0-0.09% ara', '5-10 mM IPTG, 0-0.09% ara', ...
    '0-1 mM IPTG, 0.1-0.19% ara', '1-4.9 mM IPTG, 0.1-0.19% ara',...
    '5-10 mM IPTG, 0.1-0.19% ara', '0-1 mM IPTG, 0.2% ara', ...
    '1-4.9 mM IPTG, 0.2% ara', '5-10 mM IPTG, 0.2% ara'});

new_classlabels = cell(1, length(classLabels));
for i = 1:length(new_classlabels)
    new_classlabels{i} = string(classLabels(i));
end

figure(2);
clf('reset');
% plotconfusion(y, est_cats);
% cm = confusionchart(y_ord, est_cats);
cm = confusionchart(double(comboGeneAUCData.y_ord), double(comboGeneAUCData.est_cats));
% titletext = strcat("Combo Sensor Prediction Accuracy binning into ", ...
%     strjoin(string(classLabels), ", "));
titletext = {"Combo Sensor Prediction Accuracy binning into: ", ...
    new_classlabels{:}};
% cm.GridVisible = 'off';
% Note, we can't change the grid color in here. would have to use
% illustrator
title({titletext{:}, sprintf("mAUC: %g", comboGeneAUCData.max_auc)});
ax = gca;
formatLinePlot(ax);
cm.FontSize = 12;

% Save the figure
set(gcf, 'Color', 'none');
set(gcf, 'WindowStyle', 'normal'); % in case it was docked
savefigname = "Fig4f_CM.svg";
export_fig(gcf, fullfile(savefigdir, savefigname), '-nocrop', '-Painters');
% exportgraphics(h, fullfile(savefigdir, 'Fig4e_CVs_Heatmap.eps' ));
pause(0.1);
close(gcf);

%% Make AUC  msmts latex-style table

% Need to make a table which has the measurement used and for each, the
% weight

input_tbl = struct('data', {}, 'tableRowLabels', {}, 'tableColLabels', {});
input_tbl(1).data = [comboGeneAUCData.weights(:)];
input_tbl.tableRowLabels = comboGeneAUCData.varnames_used;

input_tbl.tableColLabels = {'Weights'};
input_tbl.dataNanString = '-';
input_tbl.tableColumnAlignment = 'c';
input_tbl.tableBorders = 1;
input.booktabs = 1;
input.makeCompleteLatexDocument = 1;
latex = latexTable(input_tbl);

% here are the fields in the table to fill in:
% input.data = [1.12345 2.12345 3.12345; ...
%               4.12345 5.12345 6.12345; ...
%               7.12345 NaN 9.12345; ...
%               10.12345 11.12345 12.12345];

% input.tableColLabels = {'col1','col2','col3'};
% input.tableRowLabels = {'row1','row2','','row4'};


% % Define how NaN values in input.tableData should be printed in the LaTex table:
% input.dataNanString = '-';
%
% % Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
% input.tableColumnAlignment = 'c';
%
% % Switch table borders on/off:
% input.tableBorders = 1;
%
% % Switch table booktabs on/off:
% input.booktabs = 1;

% % Switch to generate a complete LaTex document or just a table:
% input.makeCompleteLatexDocument = 1;
%
% % % Now call the function to generate LaTex code:
% latex = latexTable(input);

%% Figure 4: flgM Unet Plot 1: Train/Val Loss

% Import data from spreadsheet

opts = spreadsheetImportOptions("NumVariables", 13);

% Specify sheet and range
opts.Sheet = "vgg11_UNet_flgm_110721_TrainVal";
opts.DataRange = "A2:M81";

% Specify column names and types
opts.VariableNames = ["Epoch", "Train Loss", "Val Loss", "Train Accuracy", "Val Accuracy", ...
    "Train Precision", "Val Precision", "Train Recall", "Val Recall", "Train IoU", ...
    "Val IoU", "Train Fscore", "Val Fscore"];
tempvartypes = strings(1, length(opts.VariableNames));
tempvartypes(:) = "double";
opts.VariableTypes = tempvartypes(:);

% Import the data
data_file = "/Users/anjali/Downloads/vgg11_UNet_flgm_110721_TrainValcsv.xlsx";
flgm_unet_noaug_trainval = readtable(data_file, opts, "UseExcel", false);

% Convert to output type
flgm_unet_noaug_trainval_array = table2array(flgm_unet_noaug_trainval);

% Clear temporary variables
clear opts tempvartypes

% Plot the Train & Val Loss
x_epoch = flgm_unet_noaug_trainval_array(:, 1);
y_trainacc = flgm_unet_noaug_trainval_array(:, 4);
y_valacc = flgm_unet_noaug_trainval_array(:, 5);

figure(2);
clf('reset');
hold on;
plot(x_epoch, y_trainacc, 'LineWidth', 3, 'DisplayName', 'Train Accuracy');
plot(x_epoch, y_valacc, 'LineWidth', 3, 'DisplayName', 'Val Accuracy');
l = legend('Location', 'southeast');

% Now let's format
% for the colors, check w tal?

ax = gca;
% set(ax, 'XScale', 'log');
ylim([0.4 1.0]);
xlim([0 60]);
title('UNet Train & Val Accuracy')

ax.XTick = [0 20 40 60];
% ax.XTickLabel = {'0', '0.5', '0.7', '1'};
ax.YTick = [0.4 0.6 0.8 1.0];
ax.XLabel.String = 'Epoch';
ax.YLabel.String = 'Accuracy';
keeplinewidths = true;
axfontsize = 24;
smallticks = false;
titlesize = 24;
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...
    keeplinewidths, axfontsize, smallticks, titlesize);

fig = gcf;

% Things to do for saving
set(fig, 'WindowStyle', 'normal'); % in case it was docked
set(fig, 'Color', 'none');
% set(gcf, 'Position', [680   385   560   420]);
set(gcf, 'Position', [680   301   453   497]); % on laptop for narrow fig
drawnow;
export_fig(fig, fullfile(savefigdir, 'Fig4_i_TrainValAcc.svg'), '-nocrop', '-Painters');
set(fig, 'WindowStyle', 'docked');
fig.Color = [1 1 1];

%% Fig 4 i: Plot Recall

y_trainrecall = flgm_unet_noaug_trainval_array(:, 8);
y_valrecall = flgm_unet_noaug_trainval_array(:, 9);

figure(2);
clf('reset');
hold on;
plot(x_epoch, y_trainrecall, 'LineWidth', 3, 'DisplayName', 'Train Recall');
plot(x_epoch, y_valrecall, 'LineWidth', 3, 'DisplayName', 'Val Recall');
l = legend('Location', 'southeast');
% Now let's format
% for the colors, check w tal?

ax = gca;
% set(ax, 'XScale', 'log');
ylim([0.4 1.1]);
xlim([0 60]);
title('UNet Train & Val Recall')
ax.YLim = [0 1.1];
ax.XTick = [0 20 40 60];
% ax.XTickLabel = {'0', '0.5', '0.7', '1'};
% ax.YTick = [0.4 0.6 0.8 1.0];
ax.YTick = [0 0.5 1.0];
ax.XLabel.String = 'Epoch';
ax.YLabel.String = 'Accuracy';
keeplinewidths = true;
axfontsize = 24;
smallticks = false;
titlesize = 24;
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...
    keeplinewidths, axfontsize, smallticks, titlesize);

fig = gcf;

% Things to do for saving
set(fig, 'WindowStyle', 'normal'); % in case it was docked
set(fig, 'Color', 'none');
% set(gcf, 'Position', [680   385   560   420]);
set(gcf, 'Position', [680   301   453   497]); % on laptop for narrow fig
drawnow;
export_fig(fig, fullfile(savefigdir, 'Fig4_i_TrainValRecall.svg'), '-nocrop', '-Painters');
export_fig(fig, fullfile(savefigdir, 'Fig4_i_TrainValRecall.png'), '-nocrop', '-Painters');
set(fig, 'WindowStyle', 'docked');
fig.Color = [1 1 1];



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions

% TIMELAPSE PLOTTING FUNCTION
function vals = getMsmtsToPlot(allExptData, field_name, iptgs, gene_list, phase_subset, skip_expts)
    vals = cell(length(gene_list), length(iptgs));
    
    % Use phase_subset as 'all', 'middle', 'end', or 'first' to choose
    % which values to include

    % Iterate over allexptdata, add constimes to appropriate place in data
    % structure
    
    % Skip any desired expts
    for i = 1:length(allExptData)
        if any(skip_expts==i)
            continue;
        end
    % Iterate over plates
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
                    tempvals = vals{geneind, iptgind};
                    if iscell(allExptData(i).(field_name)(1))
                        tempval = allExptData(i).(field_name){j};
                        if size(tempval, 1) ~= 1
                            tempval = tempval';
                        end
                        
                    else
                        tempval = allExptData(i).(field_name)(j);
                        if size(tempval, 1) ~= 1
                            tempval = tempval';
                        end
                        
                    end
                    
                    if strcmpi(phase_subset, 'all')
                        tempval = tempval;
                    elseif strcmpi(phase_subset, 'last')
                        tempval = tempval(end);
                    elseif strcmpi(phase_subset, 'first')
                        tempval = tempval(1);
                    elseif strcmpi(phase_subset, 'middle')
                        if length(tempval)>2
                            tempval = tempval(2:(length(tempval)-1));
                        elseif length(tempval)==2
                            tempval = tempval(1);
                        else
                            tempval = tempval;
                        end
                    end
                    
                    tempvals = [tempvals, tempval];
                    vals{geneind, iptgind} = tempvals;
            end
        end
    end
end

function [mag_ifftim, phase_ifftim, real_ifftim, imag_ifftim] = getIFFT(newim)
    % compute inverse fourier transform
    ifftim = ifft2(newim);
    ifftim = ifftshift(ifftim);
    mag_ifftim = abs(ifftim);
    phase_ifftim = angle(ifftim);
    real_ifftim = real(ifftim);
    imag_ifftim = imag(ifftim);
end

function [mag_fftim, phase_fftim, real_fftim, imag_fftim] = getFFT(newim)
    fftim = fft2(newim);
    fftim = fftshift(fftim);
%     mag_fftim = log(1+abs(shift_fftim));
    mag_fftim = abs(fftim);
    phase_fftim = angle(fftim);
    real_fftim = real(fftim);
    imag_fftim = imag(fftim);
end

function formatLinePlot(ax_all, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
%     ax.FontName = 'Myriad';

    if ~exist('plotgrid', 'var')
        plotgrid = false;
    end
    
    
    for axind = 1:length(ax_all)
        ax = ax_all(axind);
        if ~strcmpi(get(ax_all, 'Type'), 'ConfusionMatrixChart') &...
                ~strcmpi(get(ax_all, 'Type'), 'heatmap')

            ax.Box = 'off';
            ax.LineWidth = 2;
            ax.XColor = [0.2 0.2 0.2]; 
            ax.YColor = [0.2 0.2 0.2];
            ax.TickDir = tickdirval;    
            ax.TickLength = ticklengthval;

            if smallticks
                ax.FontSize = 14;
                if ~isempty(ax.XLabel)
                    ax.XLabel.FontSize = 24;
                end
                if ~isempty(ax.YLabel)
                    ax.YLabel.FontSize = 24;
                end
            else
                % Change both x & y labels, and tick mark fonts
                ax.FontSize = axfontsize;
            end

            % get lines
            if ~keeplinewidths
                if length(ax.Children) == 1
                    linetype = ax.Children.Type;
                    if strcmpi(linetype, 'errorbar')
                        errorlines = ax.Children;
                        errorlines.LineWidth = 3;
                        errorlines.CapSize = 1;
                    elseif strcmpi(linetype, 'line')
                        lines = ax.Children;
                        lines.LineWidth = 3;
                    end
                elseif length(ax.Children)> 1
                    for j = 1:length(ax.Children)
                        temp_part = ax.Children(j);
                        temptype = ax.Children(j).Type;
                        if strcmpi(temptype, 'errorbar')
                            temp_part.LineWidth = 3;
                            temp_part.CapSize = 1;
                        elseif strcmpi(temptype, 'line')
                            temp_part.LineWidth = 3;
                        elseif strcmpi(temptype, 'scatter')
                            % for now don't do anything?

                        end
                    end
                end
            end


            % Check if legend, if so, format
            if ~isempty(ax.Legend)
                ax.Legend.Color = [1 1 1];
                ax.Legend.EdgeColor = [0.9 0.9 0.9];
                ax.Legend.FontSize = 12;
                ax.Legend.LineWidth = 2;
                ax.Legend.TextColor = ax.XColor;
            end

            if ~isempty(ax.Title)
                ax.Title.Color = ax.XColor;
                ax.Title.FontSize = titlesize;
                ax.Title.FontWeight = 'bold';
            end

            % These options can be changed if we want a background color/grid
            if plotgrid
                ax.Color = [0.9098    0.9098    0.9098];
                ax.XGrid = 'on'; 
                ax.YGrid = 'on';
                ax.GridColor = [0.95 0.95 0.95];
                ax.GridAlpha = 1;
                ax.MinorGridColor = [0.95 0.95 0.95];
                ax.MinorGridLineStyle = '-';
                ax.MinorGridAlpha = 1;
            end
        end
    end
end

function formatHeatmapTrajPlot(ax, tickdirval, ticklengthval, heatmap_to_use, ...
    axesvis, rulersvis, titlesize)
    
    ax.Box = 'off';
    drawnow;
    pause(0.1);
    
    if axesvis & rulersvis
        ax.LineWidth = 2;
        ax.XColor = [0.2 0.2 0.2]; 
        ax.YColor = [0.2 0.2 0.2];
        ax.TickDir = tickdirval;    
        ax.TickLength = ticklengthval;
    elseif axesvis & ~rulersvis
        ax.YRuler.Axle.Visible = 'off';
        ax.XRuler.Axle.Visible = 'off';
        ax.XColor = [0.2 0.2 0.2]; 
        ax.YColor = [0.2 0.2 0.2];
        ax.TickDir = tickdirval;    
        ax.TickLength = ticklengthval;
    else
        ax.XAxis.Visible = 'off';
        ax.YAxis.Visible = 'off';
    end
    
    
    if ~isempty(ax.XLabel)
        ax.XLabel.FontSize = 24;
    end
    if ~isempty(ax.YLabel)
        ax.YLabel.FontSize = 24;
    end
    if ~isempty(ax.Title)
        ax.Title.Color = ax.XColor;
        ax.Title.FontSize = titlesize;
        ax.Title.FontWeight = 'bold';
    end
    ax.FontSize = 24;

end


function formatHeatMapChart(h, h_colormap_to_use)
    h.FontSize = 24;
    h.GridVisible = 'off';
    h.Colormap = h_colormap_to_use;
    h.CellLabelFormat = '%0.2g';
%     h.CellLabelColor = [0 0 0.1];
end


function formatBarPlot(h, ax, ticklength)
    % h is the bar plot
    if length(h) == 1
        h.EdgeColor = 'none';
    else
        for i = 1:length(h)
            h(i).EdgeColor = 'none'; 
        end
    end
    % ax is the axes
    ax.FontSize = 24;
    ax.Box = 'off';
    ax.LineWidth = 2;
    ax.XColor = [0.2 0.2 0.2]; 
    ax.YColor = [0.2 0.2 0.2];
    
    if ~isempty(ax.XLabel)
        ax.XLabel.FontSize = 24;
    end
    if ~isempty(ax.YLabel)
        ax.YLabel.FontSize = 24;
    end
    if ~isempty(ax.Title)
        ax.Title.Color = ax.XColor;
        ax.Title.FontSize = 24;
        ax.Title.FontWeight = 'bold';
    end
    
    ax.TickLength = ticklength;
    
    drawnow;
    pause(0.1);

end

function formatConfusMatrix(cm, cm_heatmap_to_use)

    cm.FontSize = 14;
    cm.DiagonalColor = cm_heatmap_to_use(end, :);
    cm.OffDiagonalColor = cm_heatmap_to_use(1, :);
    cm.GridVisible = 'off';


end

function formatLegendText(l, gene_list, linecolors)
    curr_string = l.String;
    for i = 1:length(curr_string)
        tempstring = curr_string{i};
        tempgeneind = find(strcmpi(tempstring, gene_list));
        tempcolor = linecolors(tempgeneind, :);
        newstr = sprintf('\\color[rgb]{%f, %f, %f}%s', tempcolor(1), tempcolor(2), ...
            tempcolor(3), tempstring);
        l.String{i} = newstr;
    end
    l.Interpreter = 'tex';

end

function tempvals = mkExptCell(allExptData, i, fnames)
    numims = length(allExptData(i).imfiles);
    tempvals = cell(numims, length(fnames));
    
    for j = 1:length(fnames)
        tempdata = allExptData(i).(fnames{j});
        % Get the length
        if ~isstr(tempdata) & length(tempdata)>1
        % If the length is > 1, iterate over to populate the cell
            for k = 1:(numims)
                if ~iscell(tempdata)
                    tempvals{k, j} = tempdata(k);
                else
                    tempvals{k, j} = tempdata{k};
                end

            end
        end
    end
end

