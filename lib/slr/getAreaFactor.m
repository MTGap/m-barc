function [eta] = getAreaFactor(n, thetaInc)
    % Calculate the effective area reduction factor
    % for a given incidence angle
    %
    % n: Index of refraction
    % thetaInc: Incidence angle (radians)
    %
    % getAreaFactor(1.455, pi/10) 
    %
    % Michael Gapczynski
    %------------------------------------------------
    
    % Remove nonphysical incidence angles
    thetaInc(thetaInc > pi/2) = NaN;

    thetaRef = asin(sin(thetaInc)./n);      % Internal refracted angle   
    mu = sqrt(1-tan(thetaRef).^2);
    eta = (2/pi)*(asin(mu)-sqrt(2)*tan(thetaRef)).*cos(thetaInc);
    eta(isnan(eta)) = 0;

    % Remove nonphysical reduction factors
    eta(eta < 0) = 0;
end
