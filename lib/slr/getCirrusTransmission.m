function [Tc] = getCirrusTransmission(thetaZen, t)
    % Calculate the one-way cirrus transmittance at a given zenith angle
    %
    % thetaZen: Zenith angle (radians)
    % t: Mean cirrus cloud thickness (km) - Optional
    %
    % cirrusTransmission(pi/6)
    %
    % [Degnan, Millimeter Accuracy Satellite Laser Ranging: A Review]
    %
    % Michael Gapczynski
    %--------------------------------------------------------------------
    
    if nargin < 2
        t = 1.341;      % Mean cirrus cloud thickness (km)
    end
    
    Tc = exp(-0.14.*(t.*sec(thetaZen)).^2);
    
    % "Sub-visible cirrus clouds are overhead about 50% of the time at most
    % locations. A gload study of cirrus cloud thickness [Hall et al., 1983]
    % yields the histogram in Figure 9a. From the histogram, one can compute a
    % mean cirrus cloud thickness, when present, of 1.341 km." - Degnan
end
