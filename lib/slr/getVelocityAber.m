function [alpha] = getVelocityAber(p, thetaZ, w)
    % Calculate the velocity aberration
    % 
    % p: Distance to s/c from geocenter
    % thetaZ: Zenith angle (radians)
    % w: acos((p X r)v)
    %
    % Michael Gapczynski
    %----------------------------------------------
    
    c = 299792458;    % Speed of light (m/s)
    g = 9.80665;      % Acceleration due to gravity (m/s^2)
    Re = 6371e3;      % Mean radius of Earth (m)
    
    alphaMax = (2/c)*sqrt((g*Re^2)./(p));
    alpha = alphaMax.*sqrt(cos(w).^2 + (1-((Re.*sin(thetaZ))./(p)).^2).*sin(w).^2);
end
