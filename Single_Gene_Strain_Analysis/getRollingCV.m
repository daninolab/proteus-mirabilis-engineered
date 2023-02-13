function [cvs_all, cv_traj] = getRollingCV(polarim, winwidth)
    % This function, given a window width to use, will get windows of that
    % sliding from left to right for each row of the image and get the CV
    % of each window (ie, a local CV).
    polarim = im2double(polarim);
    pol_size = size(polarim);
    cvs_all = zeros(pol_size(1), pol_size(2)-winwidth+1);
    
    % Iterate
    for col_num = (winwidth+1):pol_size(2)
        imwin = polarim(:, (col_num-winwidth):col_num);
        cvs_all(:, col_num-winwidth) = std(imwin, 0, 2)./(mean(imwin, 2));
    end
    cv_traj = mean(cvs_all, 2);
    


end