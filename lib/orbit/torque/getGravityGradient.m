function [tau] = getGravityGradient(sc, t, dt, q, w, e, nu, nudot)
    % Calculate the torque exerted on the s/c from gravity gradient
    %
    % sc: s/c configuration
    % t: Time vector (s)
    % dt: Time step (s)
    % q: Quaternion
    % w: Angular velocities in body frame
    % e: Eccentricity
    % nu: True anomaly (radians)
    % nudot: Rate of change of true anomaly (radians/s)
    %
    % Michael Gapczynski
    %--------------------------------------------------------------

    % s/c moments of inertia
    J1 = sc.J(1);
    J2 = sc.J(2);
    J3 = sc.J(3);

    if isscalar(t)
        % Determine index to look at for time
        i = ceil(t/dt);
        if i == 0
            i = 1;
        end
        e = e(i);
        nu = nu(i);
        nudot = nudot(i);
    end

    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);

    DCM = [1-2*(q2^2+q3^2) 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4);
            2*(q1*q2-q3*q4) 1-2*(q1^2+q3^2) 2*(q2*q3+q1*q4);
            2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) 1-2*(q1^2+q2^2)];

    tau = ((3*nudot^2)/(1+e*cos(nu)))*[(J3-J2)*DCM(2,1)*DCM(3,1);
                                      (J1-J3)*DCM(1,1)*DCM(3,1);
                                      (J2-J1)*DCM(1,1)*DCM(2,1)];
end
