%% calc_boundary_slopes

function [std_boundary_locations] = calc_boundary_slopes(locations, num_boundaries_mode)

    std_boundary_locations = zeros(num_boundaries_mode,1);
    
    for row = 1:num_boundaries_mode
        this_row = locations(row,:);
        this_row_reduced = this_row(find(this_row)); % ignore the 0s
        row_std = std(this_row_reduced);
        std_boundary_locations(row,1) = row_std;
    end

end