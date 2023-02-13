%% Plot Marian flgm temp sensitivity data
load("/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/MS_flgM_temperature_sensitivity/old/flgM_TempSensitivity_mean_std_sem.mat");

% plot
figure(1);
clf('reset');
hold on;
tempcolors = linspecer(4);
tempcolors = tempcolors(3:4, :);

temp_vals = [34, 36, 37];
temperature_vals = {"37C", "36C", "34C"};
iptg_vals = {"_0i_", "_10i_"};
x_vals = [0.9, 1.1, 2.9, 3.1, 4.9, 5.1];
data_names = [new_struct.Original_Msmts_File];
x_count = 1;
for i = 1:3
    temp_val = temperature_vals{i};
    disp(temp_val);
    for j = 1:2
        temp_iptg = iptg_vals{j};
        disp(temp_iptg);
        temp_ind = find(contains(data_names, temp_val) &...
            contains(data_names, temp_iptg));
%         plot(x_vals(i*j), new_struct.OverallMeanWidth_CM(temp_ind));
        scatter_vals = new_struct(temp_ind).MeanWidths_CM{1};
        scatter_x = repelem(x_vals(x_count), length(scatter_vals));
        % show the N of plates
        disp(length(scatter_vals));
        if contains(temp_iptg, "10")
            scatter(scatter_x, scatter_vals, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', tempcolors(1, :));
        else
            scatter(scatter_x, scatter_vals, 'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', tempcolors(2, :));
        end
        
        errorbar(x_vals(x_count), new_struct(temp_ind).OverallMeanWidth_CM,...
            new_struct(temp_ind).OverallSEMWidth_CM, ...
            'Marker', 'square') %, ...
%             'MarkerFaceColor','auto',...
%             'MarkerEdgeColor', 'none', 'MarkerSize', 50);
        x_count = x_count + 1;
    end
end
title("IPTG readout at different temperatures");


ax = gca;
% formatLinePlot(ax);
ax.XLim = [0 6];
ax.XTick = [1 3 5];
ax.XTickLabel = {'37', '36', '34'};
ax.YLim = [0.2 0.65];
ax.XLabel.String = 'Temperature (C)';
ax.YLabel.String = 'Mean Width of First Ring (cm)';
% ax.YTick = [1 2 3 4];
ticklengthval = [0 0];
tickdirval = 'out';


keeplinewidths=true;
axfontsize = 24; smallticks = false; 
titlesize = 24;
plotgrid = false;
legend('', '0 IPTG', '', '10 IPTG', '', '0 IPTG', '', '10 IPTG', '', '0 IPTG', '',...
    '10 IPTG', ...
    'Location', 'best');

all_e = ax.Children;

for i = 1:length(all_e)
    temp = all_e(i);
    if strcmpi(temp.Type,'errorbar') 
%         disp(temp.DisplayName);     
        temp.LineWidth = 2;
        temp.CapSize = 20;
        temp.MarkerSize = 10;
        if contains(temp.DisplayName, "10")
            % 10 IPTG color
            temp.MarkerFaceColor = tempcolors(1, :);
            temp.MarkerEdgeColor = tempcolors(1, :);
            temp.Color = tempcolors(1, :);
        else
            % 1 IPTG color
            temp.MarkerFaceColor = tempcolors(2, :);
            temp.MarkerEdgeColor = tempcolors(2, :);
            temp.Color = tempcolors(2, :);
        end

        
    end

end
ax.FontName = 'Myriad Pro';
formatLinePlot(ax, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
%% Save the figure
fig = gcf;
savefigdir = '/Users/anjali/Dropbox/_Shared_Patterning/Revision_Draft_Figures/Matlab_Plots/';
% Things to do for saving
set(fig, 'WindowStyle', 'normal'); % in case it was docked
set(fig, 'Color', 'none');
set(gcf, 'Position', [680   385   560   420]);
% set(gcf, 'Position', [217   206   708   549]); % on laptop for narrow fig
drawnow;
export_fig(fig, fullfile(savefigdir, 'FigS14_ringwidths.svg'), '-nocrop', '-Painters');

fig.Color = [1 1 1];
set(fig, 'WindowStyle', 'docked');

%% Formatting function

function formatLinePlot(ax_all, plotgrid, tickdirval, ticklengthval, ...     
    keeplinewidths, axfontsize, smallticks, titlesize);
    ax.FontName = 'Myriad Pro';

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
                        errorlines.CapSize = 5;
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