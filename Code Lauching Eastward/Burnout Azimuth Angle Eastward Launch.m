%% Determine (i) the burnout azimuth of the spacecraft for the Technical Case East Launch
%Starting Location: Kennedy's Space Station in cape canaveral
%Final Destination: coastline of california (Prep for re-entry)



clear
clc


%defining parameters

a = 6730.525196784;                  %semi-major axis (km)
e = 0.0219118;                  %Orbit's eccentricity
theta1 = 23.969*pi/180;         %initial true anomally at burnout
GSMT = 0;                       %Greenwich Siderial meridian time at burnout
mu_earth = 398600.4415;         %Gravitational parameter of earth (km^3/s^2)
R_earth = 6378.1363;            %Radius of earth (km)
J2 = 0.0010826267;              %J2 perturbations due to earths's oblaqueness
w_earth = 2*pi / 86164;         %Earth's rotation rate
n = 3;
iterations = 3;
i = 1;                         %counting variable
f = 3.353*10^-3;               %Earth's Oblaqueness Constant
height = 0;


% defining initial and final longitude/lattitude

sat_longitude1 = 279.45*pi/180;          %longitude at burnout
sat_lattitude1 = 28.50*pi/180;           %lattitude at burnout

sat_longitude_final = 241*pi/180;      %final longitude over desire position
sat_lattitude_final = 34*pi/180;       %final lattitude over desire position

% intermediate calculations
P = 2*pi*sqrt(a^3/mu_earth);
y = sqrt(1-e)*sin(theta1/2);
x = sqrt(1+e)*cos(theta1/2);
E0 = 2*atan2(y,x);
t1 = (P/(2*pi))*(E0 - e*sin(E0));      %initial time at burnout true anomally


% initiate iterations
lamnda2e = sat_longitude_final + n*w_earth*P;
deltaT = (P/(2*pi)) * (lamnda2e - sat_longitude1);   %change in time for 1st iteration only


while i <=  iterations
 
    delta_longitude = ( lamnda2e - sat_longitude1 + w_earth*deltaT );
    argument = sin(sat_lattitude_final)*sin(sat_lattitude1) + cos(sat_lattitude_final)*cos(sat_lattitude1)*cos(delta_longitude);
    theta2 = ((acos(argument) ) + theta1 );     %computes the true anomally at the next closest true anomally value to final destination

    %updating delta t
   y2 = sqrt(1-e)*sin(theta2/2);
   x2 = sqrt(1+e)*cos(theta2/2);
   E2 = 2*atan2(y2,x2);
   t2 = (P/(2*pi))*(E2 - e*sin(E2));            %new time at theta2
   
    
   

    %Computing Azimuth 
    delta_theta = theta2 - theta1;      %computes the change in anomally based on final minus initial true anomally
    argument2 = sin(delta_longitude)*cos(sat_lattitude_final) / sin(delta_theta);
    azimuth = asin(argument2);
    inclination = acos(sin(azimuth)*cos(sat_lattitude1));
   
      %accounting for Earth's oblaqueness
     if i == 1
        argument_periapsis=atan2(tan(sat_lattitude1),abs(cos(azimuth)))-theta1;
    
 
    
        constant1 = sqrt(mu_earth);
        constant2 = (1- (e^2 ))^2;
        constant3 = (a)^(7/2);
        constant4 = (sin(inclination) )^2;
        RAAN_dot = (-3*((sqrt(mu_earth)*J2*(R_earth)^2)*cos(inclination)))/(2*((1-e^2)^2)*(a^(7/2)));
        W_dot = (-3/2)*(sqrt(mu_earth)*J2*(R_earth^2)/( (1-(e^2))^2 *(a^(7/2))))*((5/2)*((sin(inclination))^2)-2);
        
        
        delta_RAAN = RAAN_dot*(deltaT   + n*P);  %fix delta_RAAN it is giving you a rounding error somewhere
        delta_W = W_dot*(deltaT + n*P);        %fix delta_W it is giving you a rounding error somewhere
        
        delta_lattitudePM = ((sin(inclination)*cos(argument_periapsis + theta2)/(cos(sat_lattitude_final)))*delta_W);
    
        %intermediate calculations to find delta_LongitudePM
        denominator = (cos(argument_periapsis + theta2) )^2;
        secant_constant = 1/denominator;
       
        num=cos(inclination)*(sec(argument_periapsis + theta2))^2;
        den=1+((cos(inclination))^2)*((tan(argument_periapsis + theta2))^2);
        delta_longitudePM=((num/den)*delta_W) + delta_RAAN;
        %update the satellite's final longitude and geodetic lattitude
  
        sat_longitude_final = (sat_longitude_final - delta_longitudePM );
        sat_lattitude_final = (sat_lattitude_final - delta_lattitudePM);
        lamnda2e = sat_longitude_final + n*w_earth*P;
        fprintf ('Delta Lattitude %1.0f  = %1.14f degrees \n',2, delta_lattitudePM*180/pi)
        fprintf ('Delta Longitude %1.0f  = %1.14f degrees \n',2, delta_longitudePM*180/pi)
     
     end
     deltaT = t2 - t1;
     fprintf ('Azimuth %1.0f iteration at burnout = %1.13f degrees\n', i, azimuth*180/pi)
     fprintf ('Inclination %1.0f iteration at burnout = %1.13f degrees\n', i, inclination*180/pi)
     fprintf ('deltaT  %1.0f iteration at burnout = %1.14f (s) \n', i, deltaT)
      disp(' ')

      
      
     

    i = i+1;        %update counting variable 
  
   
end

%Burnout Azimuth
 fprintf ('Azimuth at burnout = %1.13f degrees\n', azimuth*180/pi)
 fprintf ('deltaT at burnout = %1.14f (s) \n',  deltaT)
 disp(' ')
%computing the orbital elements at burnout

longitudeN = atan2(sin(sat_lattitude1)*sin(azimuth),cos(azimuth) );
longitude_RAAN = sat_longitude1 - longitudeN;
RAAN = GSMT + longitude_RAAN;
w =  asin(sin(sat_lattitude1)/sin(inclination) ) - theta1 ;
orbital_elements = [a e theta1 (inclination) (w) (RAAN)];
disp('Orbital Elements at burnout: ')
fprintf ('a = %1.13f (km) \n',orbital_elements(1))
fprintf ('e = %1.13f  \n',orbital_elements(2))
fprintf ('True Anomally = %1.13f degrees \n',orbital_elements(3)*180/pi)
 fprintf ('Inclination = %1.13f degrees\n',orbital_elements(4)*180/pi)
fprintf ('Argument of Periapsis = %1.13f degrees \n',orbital_elements(5)*180/pi)
fprintf ('RAAN = %1.13f degrees \n',orbital_elements(6)*180/pi)

% computing the position and velocity using orbital elements at burnout
y = CartesianCoordinateConversion (orbital_elements(1), orbital_elements(2), orbital_elements(3),orbital_elements(4),orbital_elements(5),orbital_elements(6), mu_earth);
sat_position_burnout = y(1:3) ;  %satellites position at burnout in  (km)
sat_velocity_burnout = y(4:6);     %satellites velocity at burnout in  (km/s)
 
 % using ODE45 to propagate the orbit
 
tspan = 0:10:n*P + deltaT;    %time span from 0 seconds to 3 times the orbital period
 
y_dot_function = y_dot_J2perturbations(tspan,y,mu_earth);
odeopt = odeset('reltol', 1e-12, 'abstol', 1e-20);
[t, Y] = ode45(@y_dot_J2perturbations, tspan, y, odeopt, mu_earth);
 
y_out = Y(end,:);

fprintf('Position of the spacecraft after 3 orbital periods: [%1.2f %1.2f %1.2f] (km) \n',y_out(1:3))
fprintf('Velocity of the spacecraft after 3 orbital periods: [%1.8f %1.8f %1.8f] (km/s)\n', y_out(4:6))
 sat_position_overtime = zeros(3,length(Y));
 for i = 1:length(Y)

     sat_position_overtime(:,i) = Y(i,1:3);

 end
% generating the groundtracks
 load('earth_coastline.mat')
 [sat_lat, sat_long] = GroundTracks(tspan, sat_position_overtime, GSMT, w_earth, f);


figure
grid on
plot(earth_coastline(:, 1), earth_coastline(:, 2), 'k')
% `hold' the plot so that you can add your groundtrack to it later.
hold on
axis equal
xlim([-180, 180])
ylim([-90, 90])

hold on
scatter(sat_long, sat_lat,5)
hold on
plot(241,34,'.','MarkerSize',30)

title('Satellite Groundtrack (Technical Case')
xlabel('Longitude')
ylabel('Latitude')
% % computing the elevation

%computing the elevation
range = zeros(length(sat_long),1);
azimuth2 = zeros(length(sat_long),1);
elevation = zeros(length(sat_long),1);

GMST2 = zeros(1,length(tspan));

for i = 1:size(sat_position_overtime,2)

    GMST2(i) = GSMT + w_earth*tspan(i);
    [range(i), azimuth(i), elevation(i)] = sat_to_RAE ((sat_position_overtime(1:3,i)), GMST2(i), sat_lat(end)*pi/180, sat_long(end)*pi/180, height);
  

end


figure
hold on
grid on
plot(tspan,elevation)
title('Elevation VS time (Technical Case)')
xlabel('timespan (s)')
ylabel('Elevation ')











