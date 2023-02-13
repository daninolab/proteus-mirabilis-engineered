%% remove_outlier_pixels

function [cleaned_mask, final_cols_to_consider, final_num_cols_consider] = remove_outlier_pixels(reduced_mask, num_boundaries_mode, locations, cols_to_consider)

    cleaned_mask = reduced_mask;
    final_cols_to_consider = cols_to_consider;
    
    for r = 1:num_boundaries_mode
        
        boundary_r = locations(r,cols_to_consider);
        outliers_r = isoutlier(boundary_r);
        find_outliers = find(outliers_r);
        
        if ~isempty(find_outliers)
            num_outliers = length(find_outliers);
            
            for out = 1:num_outliers
                corresponding_col = cols_to_consider(find_outliers(out));
                cleaned_mask(:,corresponding_col) = 0;
                
                col_idx_new_array = find(final_cols_to_consider == corresponding_col);
                
                % col may have been removed earlier in loop
                if ~isempty(col_idx_new_array)
                    final_cols_to_consider(col_idx_new_array) = [];
                end
            end
        end
    end
    
    final_num_cols_consider = length(final_cols_to_consider);
    
end