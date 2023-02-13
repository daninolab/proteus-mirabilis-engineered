function [polar_img, dist_vec] = flattenColonyInterp(petri_no_rim, center, dpcm, interp_val)
    % flattenColonyInterp flattens a grayscale colony image given the center of
    % the colony, the conversion from pixels to cm, and an interpolation
    % number determining the dimensions of the output image. Using triangular
    % interpolation for speed (as opposed to griddata). The distance vector is
    % used to note the distance in cm between each pixel in the output image
    % for further analysis.
    % Based on:
    % https://stackoverflow.com/questions/12924598/examples-to-convert-image-to-polar-coordinates-do-it-explicitly-want-a-slick-m

    %Get cartesian coordinates for image, given the colony center
    X0 = center(1); Y0 = center(2);
    [Y, X, z]=find(petri_no_rim);
    X=X-X0; Y=Y-Y0;
    %Convert coordinates to cm
    X = X/dpcm;
    Y = Y/dpcm;
    %Get the polar coordinates for the image
    [theta, rho] = cart2pol(X, Y);
    %Get the min & max polar coordinates
    rmin = min(rho); tmin = min(theta);
    rmax = max(rho); tmax = max(theta);

    %Create interpolant & remap image to the new grid
    F = scatteredInterpolant(theta, rho, z, 'nearest');
    [rhoi, thetai] = meshgrid(linspace(rmin, rmax, interp_val), linspace(tmin, tmax, interp_val));
    polar_img = F(thetai, rhoi);
    polar_img = polar_img'; %Rows of rings will be at the top of image

    %Get the distance vector
    dist_vec = mean(rhoi, 1);
end