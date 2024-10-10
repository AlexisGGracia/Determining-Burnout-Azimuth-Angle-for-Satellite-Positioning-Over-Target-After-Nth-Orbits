
function [sat_lat,sat_long] = GroundTracks( t,sat_pos,GMST,w, f)
    %compute the greenwich sidereal time
    
    thetaG = zeros(size(t));
    thetaG(1) = GMST;
    
    for i = 2:length(t)   
         thetaG(i) = thetaG(i-1) + w*(t(i) - t(i-1));
    end
    %rotate the position from ECI to ECF 
    %computing the geocentric lattitude and longitude
    %converting from geocentric lattitude to geodetic lattitude
    recf = zeros(3, length(sat_pos));
    geoce_Lat = zeros(size(t));
    sat_long = zeros(size(t));
    sat_lat = zeros(size(t));
    for i = 1:length(t)
        
        Q = [cos(thetaG(i)) sin(thetaG(i)) 0; -sin(thetaG(i)) cos(thetaG(i)) 0; 0 0 1];
        recf(:,i) = Q*sat_pos(:,i);
        geoce_Lat(i) = asin(recf(3,i)/norm(recf(:,i)));
        sat_long(i) = atan2(recf(2,i), recf(1,i));
        sat_lat(i) = atan2 (tan(geoce_Lat(i)), (1-f)^2);
    end
    sat_long = sat_long*180/pi;
    sat_lat = sat_lat*180/pi;
  
  
end