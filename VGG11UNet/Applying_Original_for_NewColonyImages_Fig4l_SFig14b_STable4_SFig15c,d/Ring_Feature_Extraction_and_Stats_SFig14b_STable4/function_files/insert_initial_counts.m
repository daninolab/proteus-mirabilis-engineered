%% insert_initial_counts

function [msmtStruct] = insert_initial_counts(msmtStruct, n, num_boundary_pixels, num_boundaries_mode, vertical_count)

    msmtStruct(n).TotalBoundaryPixels = num_boundary_pixels;
    msmtStruct(n).NumberOfBoundariesMode = num_boundaries_mode;
    msmtStruct(n).NumberOfVerticalBoundaryPixels = vertical_count;
    
end