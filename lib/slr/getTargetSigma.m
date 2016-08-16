function [sigma] = getTargetSigma(D, rho, n, delta, alpha, thetaInc, lambda)
    % Calculate the optical cross-section of a fused silica corner cube
    %
    % D: Corner cube diameter (m)
    % rho: Reflectivity
    % n: Index of refraction
    % delta: Dihedral angle offset (radians)
    % alpha: Velocity aberration / off-axis angle (radians)
    % thetaInc: Incidence angle (radians)
    % lambda: Laser wavelength (m)
    %
    % getTargetSigma(38.1e-3, 0.93, 1.455, 7e-6, 28e-6, pi/10, 532.1e-9) 
    %
    % Michael Gapczynski
    %-------------------------------------------------------------------
    
    % Peak, on-axis, optical cross-section
    sigmaCC = (pi^3*rho*D^4)/(4*lambda^2);
    
    % A dihedral angle offset spoils the retroreflector
    % and splits the FFDP into 2N lobes, where N is the number
    % of cube angles offset, i.e. N = 1, 2, 3
    % This FFDP calculation assumes N = 3
    if delta > 0
        % Peak cross-section of each lobe is reduced by (2N)^2
        sigmaCC = (1/36)*sigmaCC;
        gamma = (4/3)*sqrt(6)*n*delta;
        theta = alpha - gamma;
    else
        theta = alpha;
    end
    
    % Get the effective area reduction factor
    eta = getAreaFactor(n, thetaInc);
    
    % Airy function
    r = ((pi*D)/lambda).*sin(theta);
    J1 = besselj(1, r);
    sigma = (eta.^2).*sigmaCC.*((2*J1)./r).^2; 
end
