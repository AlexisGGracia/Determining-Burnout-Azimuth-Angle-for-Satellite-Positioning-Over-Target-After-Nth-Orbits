% Given the Satellites position, GMST, latitude, longitude, and height, it computes the range, azimuth and elevation 

function [range,azimuth,elevation] = Elevation_function(sat_position, GMST,geodetic_lat, longitude, height)

   % Defining constants
    format long
    r_earth = 6378.1363;        % Radius of Earth (km)
    f = 3.353 * 10^-3;          % Flattening of Earth
    constant = (sin(geodetic_lat))^2;  % Constant to simplify the expressions
    top = r_earth * ((1 - f)^2);

    A = ((r_earth / sqrt(1 - (2 * f - f^2) * constant)) + height)*cos(geodetic_lat);
    B = (top / sqrt(1 - (2 * f - f^2) * constant)) + height;
    
    % Computing the relative position of the observers in the ECF frame and
    % converting inertial position to ECF
    
    Recf = [A * cos(longitude); A * sin(longitude); B * sin(geodetic_lat)];  % Assuming B is related to geodetic latitude for the Z-axis


    Qijk2ECF = [cos(GMST) sin(GMST) 0; -sin(GMST) cos(GMST) 0; 0 0 1];
    recf = Qijk2ECF.*sat_position;       %satellites position in the ECF frame
    pecf = recf - Recf;                  %satellites relative position in the ECF frame

    % Converting the satellites relative position from ECF to ENZ
    r2_latt = [cos(-geodetic_lat) 0 -sin(-geodetic_lat); 0 1 0; sin(-geodetic_lat) 0 cos(-geodetic_lat)];
    r3_long = [cos(longitude) sin(longitude) 0; -sin(longitude) cos(longitude) 0; 0 0 1];
    Qecf2enz = [0 1 0; 0 0 1; 1 0 0]*r2_latt*r3_long;
    Penz = Qecf2enz*pecf;          %satellites relative position in the ENZ frame (needed for 2.2)

    % 2.2 Question
        
    range = norm(Penz);
    azimuth = atan2(Penz(1),Penz(2));
    elevation = asin(Penz(3)/norm(Penz));

end

