%% cols_with_boundaries_mode

function [num_boundary_pixels, labeled_mask, num_boundaries_mode, cols_to_consider, num_cols_consider] = cols_with_boundaries_mode(this_mask)

    % number of ring boundary pixels detected (including those we will omit later)
    num_boundary_pixels = length(find(this_mask));

    % add labels to ring boundaries
    [labeled_mask, num_rings] = bwlabel(logical(this_mask));
    
    [h w] = size(labeled_mask);

    num_boundaries_all_cols = zeros(1,w);

    for c = 1:w

        % grab 1 column
        col_c = labeled_mask(:,c);

        % grab the labels
        col_c_labels = unique(col_c);

        % drop 0 (not a label, just background)
        col_c_labels(col_c_labels == 0) = [];

        % count # of ring boundaries detected
        num_boundaries = length(col_c_labels);

        num_boundaries_all_cols(1,c) = num_boundaries;
    end

    num_boundaries_mode = mode(num_boundaries_all_cols);
    
    % Prior knowledge: there are at least 2 boundaries 
    % (that of the inoculum edge, and at the least ring formed beyond the inoculum)
    if num_boundaries_mode == 1
        num_boundaries_mode = 2;
    end
    
    cols_to_consider = find(num_boundaries_all_cols == num_boundaries_mode);
    num_cols_consider = length(cols_to_consider);

end