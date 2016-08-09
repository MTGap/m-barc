function [t, pos, vel, nu] = getOrbit(hA, hP, inc, w, bigOmega, m0, dt)
    % Calculate the satellite's position and velocity of the specified orbit
    % for the provided time step
    %
    % hA: Altitude of apogee (m)
    % hP: Altitude of perigee (m)
    % inc: Inclination (radians)
    % w: Argument of perigee (radians)
    % bigOmega: Right ascension of the ascending node - RAAN (radians)
    % m0: Mean anomaly at epoch (radians)
    % dt: Time step (s)
    %
    % getOrbit(35786e3, 200e3, 28.5*(pi/180), 0, 0, 0, 100)
    %
    % Code adapted from Fundamentals of Spacecraft Attitude Determination and Control
    % TRMM Attitude Deterination Example in Chapter 5
    % Written by John L. Crassidis and F. Landis Markley, 3/31/2011
    %
    % Michael Gapczynski
    %---------------------------------------------------------------------------------

    Re = 6378e3;                % Radius of Earth (m)
    mu = 3.986e14;              % Gravitational parameter of the Earth (m^3/s^2)

    rA = Re+hA;                 % Radius of apogee (m)
    rP = Re+hP;                 % Radius of perigee (m)

    a = (rA+rP)/2;              % Semi-major axis (m)
    e = (rA-rP)/(rA+rP);        % Eccentricity
    n = sqrt(mu/a^3);           % Mean motion (radians/s)
    T = (2*pi)/n;               % Orbital period (s)
    
    t = 0:dt:T;
    Lt = length(t);

    % Attitude matrix for orbit frame conversion from perifocal to inertial
    r11 = cos(bigOmega)*cos(w)-sin(bigOmega)*sin(w)*cos(inc);
    r12 = -cos(bigOmega)*sin(w)-sin(bigOmega)*cos(w)*cos(inc);
    r21 = sin(bigOmega)*cos(w)+cos(bigOmega)*sin(w)*cos(inc);
    r22 = -sin(bigOmega)*sin(w)+cos(bigOmega)*cos(w)*cos(inc);
    r31 = sin(w)*sin(inc);
    r32 = cos(w)*sin(inc);

    % Preallocate space
    pos = zeros(Lt, 3);
    vel = zeros(Lt, 3);
    nu = zeros(Lt, 1);
    
    eSqrtPlusMinus = sqrt((1+e)/(1-e));

    % Tolerance for Newton's method iterations
    eps = 1e-10;
    count = 0;
    maxIter = 100;

    for i = 1:Lt
        % Solve Kepler's equation via Netwon's method
        m = m0+n*(t(i)-t(1));   % Mean anomaly (radians)
        bigE = m;               % Eccentric anomaly (radians)
        deltaE = 10;
        while abs(deltaE) > eps
            deltaE = (m-(bigE-e*sin(bigE)))/(1-e*cos(bigE));
            bigE = bigE+deltaE;
            count = count + 1;
            if count == maxIter, disp(' Maximum Number of Iterations Achieved'), break, end
        end
        count = 0;

        % Perifocal position
        rmag = a*(1-e*cos(bigE));
        x = a*(cos(bigE)-e);
        y = a*sqrt(1-e^2)*sin(bigE);
        xdot = -sqrt(mu*a)/rmag*sin(bigE);
        ydot = sqrt(mu*a*(1-e^2))/rmag*cos(bigE);

        % Inertial position 
        pos(i,1) = r11*x+r12*y;
        pos(i,2) = r21*x+r22*y;
        pos(i,3) = r31*x+r32*y;
        vel(i,1) = r11*xdot+r12*ydot;
        vel(i,2) = r21*xdot+r22*ydot;
        vel(i,3) = r31*xdot+r32*ydot;
        
        % True anomaly (radians)
        nu(i) = 2*atan(eSqrtPlusMinus*tan(bigE/2));    
    end
    % Add 2pi if negative
    nu = mod(nu, 2*pi);
end
