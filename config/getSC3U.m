function [sc] = getSC3U()
    % Default configuration for a 3U CubeSat
    %
    % Michael Gapczynski
    %----------------------------------------
    
    % Principal moments of inertia (kg m^2)
    J = [0.0268 0.0268 0.005];
    
    % +X face
    surf1 = struct();
    % Outward normal unit vector from surface
    surf1.normal = [1 0 0];
    % Vector from the s/c center of mass to the center of pressure on surface (m)
    surf1.cp = [0.05 0 0];
    % Surface area (m^2)
    surf1.area = (340.5*100)*1e-6;
    % Reflectance factor (0-1), 0 = perfect absorption, 1 = perfect reflection
    surf1.reflectance = 0.6;
    
    % -X face
    surf2 = surf1;
    surf2.normal = [-1 0 0];
    surf2.cp = [-0.05 0 0];
    
    % +Y face
    surf3 = surf1;
    surf3.normal = [0 1 0];
    surf3.cp = [0 0.05 0];
    
    % -Y face
    surf4 = surf1;
    surf4.normal = [0 -1 0];
    surf4.cp = [0 -0.05 0];
    
    % +Z face
    surf5 = struct();
    surf5.normal = [0 0 1];
    surf5.cp = [0 0 0.17];
    surf5.area = (100*100)*1e-6;
    surf5.reflectance = 0.6;
    
    % -Z face
    surf6 = surf5;
    surf6.normal = [0 0 -1];
    surf6.cp = [0 0 -0.17];
    
    surfaces = [surf1 surf2 surf3 surf4 surf5 surf6];
    
    sc = struct('name', '3U', 'J', J, 'surfaces', surfaces);
end

