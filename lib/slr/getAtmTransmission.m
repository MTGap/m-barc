function [Ta] = getAtmTransmission(thetaZen, hGS, V)
    % Calculate the one-way atmospheric transmission for a wavelength of
    % 532 nm at a given zenith angle
    %
    % thetaZen: Zenith angle (radians)
    % hGS: Ground station altitude above sea level (m)
    % V: Sea level visibility (1-5: Extremely Clear - Light Haze) - Optional
    %
    % atmTransmission(pi/6, 600)
    %
    % [Degnan, Millimeter Accuracy Satellite Laser Ranging: A Review]
    %
    % Michael Gapczynski
    %----------------------------------------------------------------------
    
    hScale = 1.2;   % Scale height (km)
    
    % Attenuation coefficents for 532 nm at sea level (km^-1)
    sigma = [0.0439     % Extremely Clear (V = 60 km)
             0.0930     % Very Clear (V = 40 km)
             0.1486     % Standard Clear (V = 23.5 km)
             0.2491     % Clear (V = 15 km)
             0.4748];   % Light Haze (V = 8 km)
    
    if nargin < 3
        V = 3;
    end

    Ta = exp(-sigma(V)*hScale.*sec(thetaZen).*exp(-(hGS/1000)/hScale));
end