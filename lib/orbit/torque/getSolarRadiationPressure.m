function [tau] = getSolarRadiationPressure(sc, t, dt, q, w, pos, sun)
    % Calculate the torque exerted on the s/c from solar radiation pressure
    %
    % sc: s/c configuration
    % t: Time vector (s)
    % dt: Time step (s)
    % q: Quaternion
    % w: Angular velocities in body frame
    % pos: s/c position vector in ECI frame
    % sun: Sun unit vector in ECI frame
    %
    % Michael Gapczynski
    %----------------------------------------------------------------------

    AU = 1.4960e11;     % Astronomical unit (m)
    Fs = 1367;          % Solar constant (W/m^2)
    c = 299792458;      % Speed of light (m/s)

    if isscalar(t)
        % Determine index to look at for time
        i = ceil(t/dt);
        if i == 0
            i = 1;
        end
        pos = pos(i,:);
        sun = sun(i,:);
    end

    tau = zeros(3, 1);

    % Check if s/c is in sunlight
    if inSunlight(pos, AU*sun)
        for i = 1:numel(sc.surfaces)
            surf = sc.surfaces(i);
            % Check if surface is in sunlight
            if surf.normal*sun' > 0
                tau = tau + (Fs/c)*surf.area*(1+surf.reflectance)*surf.normal*sun'.*cross(sun', surf.cp');
            end
        end
    end
end
