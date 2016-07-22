function [PSolar, solarArray] = getSolarPower(solarArray, date, dt, tf, pos, q)
    % Calculate the power produced by a solar array
    %
    % solarArray: struct of solar array panels
    % date: Start time (datetime)
    % dt: Time step (s)
    % tf: Run time (s)
    % pos: s/c position
    % q: s/c attitude
    %
    % getSolarPower(getSolarArrayZ(), datetime(2019, 1, 20), 100, 10000, pos, q)
    %
    % Michael Gapczynski
    %---------------------------------------------------------------------------
    
    P0 = 1360.8;        % Solar constant (W/m^2)
    AU = 1.4960e11;     % Astronomical unit (m)

    % Sun Reference and Eclipse
    sun = getSunEphemeris(date, dt, tf);
    sunB = i2b(q, sun);
    timeInSun = inSunlight(pos, AU*sun);
    
    PSolar = zeros(length(pos), 1);
    
    for i = 1:numel(solarArray)
        sp = solarArray{i};
        spN = repmat(sp.normal, length(pos), 1);
        c = dot(sunB, spN, 2);
        sp.thetaInc = atan2(vecnorm(cross(sunB, spN, 2)), c);
        % Remove incidence angles that can't be illuminated
        sp.thetaInc(sp.thetaInc > pi/2) = NaN;
        % Remove incidence angles during eclipse
        sp.thetaInc(timeInSun == 0) = NaN;
        sp.power = sp.eta*P0*sp.area.*cos(sp.thetaInc);
        sp.power(isnan(sp.power)) = 0;
        PSolar = PSolar + sp.power;
        solarArray{i} = sp;
    end
end
