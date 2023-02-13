%% distance_btwn_boundaries

function [overall_min_max_mean_dist, means_stds_per_dist] = distance_btwn_boundaries(num_boundaries_mode, final_num_cols_consider, final_cols_to_consider, final_locations)

    % set up matrix for storing distances between consecutive ring boundaries
    distances = zeros(num_boundaries_mode, final_num_cols_consider);
    
    for i = 1:final_num_cols_consider
        
        % set up array for storing distances btwn ring boundaries, for this column
        hit_diffs = zeros(num_boundaries_mode, 1);
        
        % grab column
        col_idx = final_cols_to_consider(i);
        ring_hits_this_col = final_locations(:,col_idx);
        
        % store the first hit in first row (i.e. radius of inoculum)
        hit_diffs(1,1) = ring_hits_this_col(1);
        
        for h = 2:num_boundaries_mode
            diff_h = ring_hits_this_col(h) - ring_hits_this_col(h-1);
            hit_diffs(h,1) = diff_h;
        end

        distances(:,i) = hit_diffs;
    end
    
    overall_min_dist = min(distances,[],'all');
    overall_max_dist = max(distances,[],'all');
    overall_mean_dist = mean(distances,'all');
    overall_min_max_mean_dist = [overall_min_dist, overall_max_dist, overall_mean_dist];
    
    % avg across the cols to get avg distance between consecutive ring boundaries for the whole image
    % and the standard deviation across the columns
    means_stds_per_dist = zeros(num_boundaries_mode, 2); % col 1 = means, 2 = stds

    for a = 1:(num_boundaries_mode)
        mean_a = mean(distances(a,:),'all');
        means_stds_per_dist(a,1) = mean_a;

        std_a = std(distances(a,:));
        means_stds_per_dist(a,2) = std_a; 
    end
end