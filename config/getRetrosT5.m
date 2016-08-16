function [retros] = getRetrosT5()
    % A configuration of 5 corner cubes on different faces
    %
    % Michael Gapczynski
    %------------------------------------------------------

    name = '±X1±Y1+Z1_38S';
    
    % Fused silica corner cube
    D = 38.1e-3;           % Diameter (m)
    rho = 0.93;            % Reflectivity
    n = 1.455;             % Index of refraction
    delta = 6e-6;          % Dihedral angle offset (radians)
    
    retro = struct('D', D, 'rho', rho, 'n', n, 'delta', delta);
    
    % +Z face
    retro1 = retro;
    retro1.normal = [0 0 1];

    % +X face
    retro2 = retro;
    retro2.normal = [1 0 0];
    
    % -X face
    retro3 = retro;
    retro3.normal = [-1 0 0];
    
    % +Y face
    retro4 = retro;
    retro4.normal = [0 1 0];
    
    % -Y face
    retro5 = retro;
    retro5.normal = [0 -1 0];
    
    retros = struct('name', name, 'retros', [retro1 retro2 retro3 retro4 retro5]);
end
