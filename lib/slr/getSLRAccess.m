function [gs] = getSLRAccess(gs, retros, pos, vel, q)
    % Find possible access times for SLR
    %
    % gs: Ground station structure array
    % retros: Retroreflector configuration structure array
    % pos: s/c position vector in ECI frame
    % vel: s/c velocity vector in ECI frame
    % q: s/c attitude quaternion
    %
    % Michael Gapczynski
    %--------------------------------------------------------
    
    % Set a maximum range for SLR attempts
    allowedRange = 20000e3;

    % Vector from ground station to s/c
    r = pos - gs.pos;
    gs.range = vecnorm(r);
    % Unit vector from geocenter to ground station
    gsNormal = gs.pos./repmat(vecnorm(gs.pos), 1, 3);
    % Calculate zenith angle for each s/c position
    gs.angle = getAngleBetween(r, gsNormal);
    
    % Find time indices of SLR access times
    inRange = find(gs.range <= allowedRange);
    inSight = find(gs.angle <= gs.maxZenith);
    gs.indices = intersect(inRange, inSight);
    
    % Compute ranging values for in range/sight times only
    posNormal = vecnorm(pos(gs.indices,:));
    p = pos(gs.indices,:)./repmat(posNormal, 1, 3);
    rNormal = vecnorm(r(gs.indices,:));
    r = r(gs.indices,:)./repmat(rNormal, 1, 3);
    velNormal = vecnorm(vel(gs.indices,:));
    v = vel(gs.indices,:)./repmat(velNormal, 1, 3);
    gs.w = acos(dot(cross(p, r, 2), v, 2));
    gs.alpha = getVelocityAber(posNormal, gs.angle(gs.indices), gs.w);

    % Transform ground station to s/c vector from ECI frame to body frame
    rB = i2b(q(gs.indices,:), r);

    Li = length(gs.indices);
    gs.sigma = zeros(Li, 1);
    gs.thetaInc = zeros(Li, numel(retros.retros));
    for i = 1:numel(retros.retros)
        retro = retros.retros(i);
        retroNormal = repmat(retro.normal, Li, 1);
        % Determine angle of each retroreflector w.r.t laser
        thetaInc = getAngleBetween(rB, retroNormal);
        gs.thetaInc(:,i) = thetaInc;
        % Sum optical cross-section from each retroreflector in view
        gs.sigma = gs.sigma+getTargetSigma(retro.D, retro.rho, retro.n, retro.delta, gs.alpha, thetaInc, gs.lambda);
    end
    % Calculate transmission losses
    Ta = getAtmTransmission(gs.angle(gs.indices), gs.LLA(3));
    Tc = getCirrusTransmission(gs.angle(gs.indices));
    % Count mean number of photoelectrons detected
    gs.n = countPhotons(gs, gs.range(gs.indices), gs.sigma, Ta, Tc);
    % Probability
    gs.Pdb = 1-exp(-gs.n);
end
