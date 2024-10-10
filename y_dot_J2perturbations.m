%Function uses two-body dynamics including J2 pertubations output dy which contains the six constants such as velocity(1:3) and acceleration(4:6)
% y is a vector that contains the three components in the XYZ of position(1:3) and velocity(4:6)
function dy = y_dot_J2perturbations(t, y,mu)
    %defining constants
    format longG
    J2 = 0.00108248;        %J2 perturbations constant (unitless)
    R = 6378.1363;          %Radius of earth (km)
    
    dy = zeros(6,1);
    r = norm(y(1:3));
    dy(1:3) = y(4:6);

    constant1 = J2*(3/2)*((R/r)^2);
    constant2 = 5*(y(3)/r)^2;
    constant = 1 - constant1*(constant2 - 1);
    kconstant = 1 - constant1*(constant2 - 3);
    i = (-mu*y(1)/(r^3))*constant;  % i-component of the acceleration
    j = (-mu*y(2)/(r^3))*constant;  % j-component of the acceleration
    k = (-mu*y(3)/(r^3))*kconstant;
    dy(4:6) = [i;j;k];
end
