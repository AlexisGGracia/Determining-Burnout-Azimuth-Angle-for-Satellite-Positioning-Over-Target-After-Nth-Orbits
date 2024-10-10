% Function converts Orbital Elements (OE) to Position and Velocity (PV)
function motion = CartesianCoordinateConversion (a, e, theta, i, w, omega, mu)
    format longG
        E = 0.00000000001;      %some very small number for comparasion purpose
        %Define the matrices to build transformation operator Q
        R1 = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
        R3_F = [cos(-omega) sin(-omega) 0; -sin(-omega) cos(-omega) 0; 0 0 1];
        R3_S = [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];
        p = a*(1-e^2);          %semi-parameter
        rpqw = zeros(3,1);
        rijk = zeros(3,1);
        vpqw = zeros(3,1);
        vijk = zeros(3,1);
        if e < E && i < E   %first case scenario
            format longG
            Q = eye(3,3);
            angle = theta;      %theta equals true longitude which is an input
        elseif e < E
            format longG
            Q = R1*R3_S;   
            angle = theta;      %theta equals argument of lattitude which is an input
        elseif i < E
            format longG
            Q = R3_F;
            angle = theta;      %theta equals true anomaly which is an input
        else
            format longG
            Q = R3_F*R1*R3_S;
            angle = theta;      %theta equals true anomaly which is an input
        end
        %compute the position and velocity vector in the pqw frame
       rpqw(1) = (p/(1+e*cos(angle)))*cos(angle);
       rpqw(2) = (p/(1+e*cos(angle)))*sin(angle);
       vpqw(1) = sqrt(mu/p)*(-sin(angle));
       vpqw(2) = sqrt(mu/p)*(e + cos(angle));
       %transform the position and velocity vectors into cartesian
       %coordinates
format longG
       rijk = Q*rpqw;
       vijk = Q*vpqw;
       motion = [rijk;vijk];
       % fprintf ('r = %1.2fi + %1.2fj + %1.2fk \n',motion(1), motion(2), motion(3))
       % fprintf ('v = %1.2fi + %1.2fj + %1.2fk \n',motion(4), motion(5), motion(6))

end
