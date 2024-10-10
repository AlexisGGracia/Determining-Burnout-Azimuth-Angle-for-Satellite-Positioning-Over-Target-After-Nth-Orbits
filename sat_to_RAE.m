% Given the Satellites position, GMST, latitude, longitude, and height, it computes the range, azimuth and elevation 

function [range, azimuth_deg, elevation_deg] = sat_to_RAE(sat_position, GMST_rad, geodetic_lat_rad, longitude_rad, height_km)

    x_unit_vector = [1 0 0];
    y_unit_vector = [0 1 0];
    z_unit_vector = [0 0 1];
    Earth_radius_km = 6378.1363;
    flattening_factor = 3.353e-3;
    
    % Calculate geodetic radius in the ECF 
    R_ECF = [(((Earth_radius_km / sqrt(1 - (2 * flattening_factor - flattening_factor^2) * (sin(geodetic_lat_rad)^2))) + height_km) * cos(geodetic_lat_rad)) * cos(longitude_rad);
             (((Earth_radius_km / sqrt(1 - (2 * flattening_factor - flattening_factor^2) * (sin(geodetic_lat_rad)^2))) + height_km) * cos(geodetic_lat_rad)) * sin(longitude_rad);
             ((((Earth_radius_km * (1 - flattening_factor)^2) / sqrt(1 - (2 * flattening_factor - flattening_factor^2) * (sin(geodetic_lat_rad)^2))) + height_km) * sin(geodetic_lat_rad))];

    % Calculate the transformation matrix from ECF to ENZ 
    QECF_ENZ = [0 1 0; 0 0 1; 1 0 0] * [cos(-geodetic_lat_rad) 0 -sin(-geodetic_lat_rad);
                0 1 0; sin(-geodetic_lat_rad) 0 cos(-geodetic_lat_rad)] * [cos(longitude_rad) sin(longitude_rad) 0;
                -sin(longitude_rad) cos(longitude_rad) 0; 0 0 1];


    QIJK_ECF = [cos(GMST_rad) sin(GMST_rad) 0; -sin(GMST_rad) cos(GMST_rad) 0; 0 0 1];

    % Transform the satellite's position from IJK to ECF
    sat_position_ECF = QIJK_ECF * sat_position;

    % Calculate the range vector in ECF coordinates
    range_ECF = -1 * R_ECF + sat_position_ECF;

    % Transform the range vector from ECF to ENZ coordinates
    range_ENZ = QECF_ENZ * range_ECF;

    % Calculate range, azimuth, and elevation
    range = sqrt(dot(range_ENZ, range_ENZ));
    azimuth_deg = (atan2(dot(range_ENZ, x_unit_vector), dot(range_ENZ, y_unit_vector))) * 180 / pi;
    elevation_deg = (asin(dot(range_ENZ, z_unit_vector) / range)) * 180 / pi;

end
