%% get_boundary_locations

function [locations] = get_boundary_locations(num_boundaries_mode, num_cols_consider, cols_to_consider, reduced_mask)

    % set up matrix for storing ring boundary locations
    [h w] = size(reduced_mask);
    locations = zeros(num_boundaries_mode, w);

    for i = 1:num_cols_consider

        % grab column
        col_idx = cols_to_consider(i);
        this_col = reduced_mask(:,col_idx);

        ring_hits = find(this_col);
        locations(:,col_idx) = ring_hits;
        
    end

end