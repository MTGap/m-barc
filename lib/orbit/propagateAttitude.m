function [q, w, tau] = propagateAttitude(sc, t, dt, q0, w0, torqueFuncs)
    % Propagate the s/c attitude with kinematic and dynamics equations
    %
    % sc: s/c configuration structure array
    % t: Time vector (s)
    % dt: Time step (s)
    % q0: Initial quaternion
    % w0: Initial angular velocities in body frame
    % torqueFuncs: Cell array of function handles to calculate torques
    %   Each row specifies a function handle with any additional arguments
    %   included in another cell array in the second column
    %   e.g. {@getGravityGradient {e nu nudot}};
    %        tau(:, 1:3) = getGravityGradient(sc, t, dt, q, w, e, nu, nudot)
    %
    % Michael Gapczynski
    %-----------------------------------------------------------------------

    % s/c moments of inertia
    J1 = sc.J(1);
    J2 = sc.J(2);
    J3 = sc.J(3);

    J2MinusJ3OverJ1 = (J2-J3)/J1;
    J3MinusJ1OverJ2 = (J3-J1)/J2;
    J1MinusJ2OverJ3 = (J1-J2)/J3;

    if exist('torqueFuncs', 'var') && iscell(torqueFuncs)
        [n, ~] = size(torqueFuncs);
    else
        n = 0;
    end

    x0 = [q0(:); w0(:)];
    options = odeset('RelTol', 1e-8);
    [~, x] = ode45(@kd, t, x0, options);
    function xdot = kd(ti, xi)
        % Quaternion
        qi = [xi(1); xi(2); xi(3); xi(4)];

        % Angular velocities
        wi = [xi(5) xi(6) xi(7)];
        w1 = wi(1);
        w2 = wi(2);
        w3 = wi(3);

        % Kinematic equations
        qdot = 0.5*[0 w3 -w2 w1;
                    -w3 0 w1 w2;
                    w2 -w1 0 w3;
                    -w1 -w2 -w3 0]*qi;

        % Dynamics equations
        wdot = [J2MinusJ3OverJ1*w2*w3;
                J3MinusJ1OverJ2*w3*w1;
                J1MinusJ2OverJ3*w1*w2];

        % Include any specified torques           
        for i = 1:n
            torqueFunc = torqueFuncs{i, 1};
            if isa(torqueFunc, 'function_handle')
                extraArgs = torqueFuncs{i, 2};
                args = [{sc ti dt qi wi} extraArgs{:}];
                taui = torqueFunc(args{:});
                wdot = wdot + taui(:)./[J1; J2; J3];
            end
        end

        xdot = [qdot; wdot];
    end

    % Format outputs from ode45 solver
    % Quaternions
    q = x(:, 1:4);

    % Angular velocities
    w = x(:, 5:7);

    tau = zeros(length(t), n*3);

% TODO probably switch loops

    % Loop through torque calculations again to return as output
    for j = 1:length(t)
        for k = 1:n
            torqueFunc = torqueFuncs{k, 1};
            if isa(torqueFunc, 'function_handle')
                extraArgs = torqueFuncs{k, 2};
                args = [{sc t(j) dt q(j,:) w(j,:)} extraArgs{:}];
                tau(j, (k*3-3)+(1:3)) = torqueFunc(args{:});
            end
        end
    end
 end
