function [light] = inSunlight(pos, sun)
    % Determine if the s/c is in sunlight at a given position
    %
    % pos: s/c position vector in ECI frame
    % sun: Sun vector in ECI frame
    %
    % inSunlight([6.492e6 0.8964e6 0.4867e6], [0.2338e11 1.3557e11 0.5877e11])
    %
    % Reference:
    % Cassandra Belle VanOutryve, "A thermal analysis and design tool for
    % small spacecraft"
    %
    % Michael Gapczynski
    %--------------------------------------------------------------------------
    
    Re = 6378e3;        % Radius of Earth (m)
    
    theta1 = acos(Re./vecnorm(pos));
    theta2 = acos(Re./vecnorm(sun));
    
    % Angle between s/c position vector and Sun position vector (radians)
    psi = acos(sum(pos.*sun, 2)./(vecnorm(pos).*vecnorm(sun)));
    
    light = ones(length(pos(:,1)), 1);
    % If psi is >= thetal+theta2, the s/c is in eclipse, otherwise it's in sunlight
    i = (psi >= theta1 + theta2);
    light(i) = 0;
end

