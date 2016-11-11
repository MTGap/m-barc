function [retros] = getRetrosZ4()
    % A planar array of 4 retroreflectors on the +Z face 
    %
    % Michael Gapczynski
    %------------------------------------------------------

    name = '+Z4_38S';
    
    % Fused silica corner cube
    D = 38.1e-3;           % Diameter (m)
    rho = 0.93;            % Reflectivity
    n = 1.455;             % Index of refraction
    delta = 6e-6;          % Dihedral angle offset (radians)
    
    retro = struct('D', D, 'rho', rho, 'n', n, 'delta', delta);
    
    retro1 = retro;
    retro1.normal = [0 0 1];
    retro2 = retro1;
    retro3 = retro1;
    retro4 = retro1;
    
    retros = struct('name', name, 'retros', [retro1 retro2 retro3 retro4]);
end
