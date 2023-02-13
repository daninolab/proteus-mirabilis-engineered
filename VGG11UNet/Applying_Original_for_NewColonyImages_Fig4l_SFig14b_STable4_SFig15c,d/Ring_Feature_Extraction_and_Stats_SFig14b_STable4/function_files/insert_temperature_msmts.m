%% insert_temperature_msmts

function [msmtStruct] = insert_temperature_msmts(n, msmtStruct, final_num_cols_consider, overall_min_max_mean_dist, means_stds_per_dist, std_boundary_locations)

    msmtStruct(n).NumberOfColumnsConsidered = final_num_cols_consider;
    
    msmtStruct(n).OverallMinDistance = overall_min_max_mean_dist(1);
    msmtStruct(n).OverallMaxDistance = overall_min_max_mean_dist(2);
    msmtStruct(n).OverallMeanDistance = overall_min_max_mean_dist(3);
    
    msmtStruct(n).MeanDistancesBtwnBoundaries = means_stds_per_dist(:,1);
    msmtStruct(n).MeanDistanceToInoc = means_stds_per_dist(1,1);
    msmtStruct(n).MeanDistanceInocToBound2 = means_stds_per_dist(2,1);
    
    msmtStruct(n).STDofDistancesBtwnBoundaries = means_stds_per_dist(:,2);
    msmtStruct(n).MaxSTDofDistances = max(means_stds_per_dist(:,2));
    msmtStruct(n).MinSTDofDistances = min(means_stds_per_dist(:,2));

    msmtStruct(n).STDof2ndBoundaryLocations = std_boundary_locations(2,1);
    msmtStruct(n).STDofLastBoundaryLocations = std_boundary_locations(length(std_boundary_locations),1);
    msmtStruct(n).MaxSTDofAllBoundaryLocations = max(std_boundary_locations); 

end