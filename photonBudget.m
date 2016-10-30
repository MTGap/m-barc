close all;

addpath(genpath('config'));
addpath(genpath('lib'));

% Photon Budget
%
% Michael Gapczynski
%----------------------------------------

sc = getSC3U();             % s/c configuration
retros = getRetrosT5();     % Retroreflector configuration

% Folder to save figures
outFolder = ['output/slr/' sc.name '/' retros.name];

scName = sc.name;
scRetrosName = [sc.name ' ' retros.name];

% Time
dt = 10;                        % Time step (s)
date = datetime(2018, 8, 20);   % Year, month, day
nOrbits = 2;                    % Number of orbits

% Orbit Parameters
hA = 35400e3;               % Altitude of apogee (m)
hP = 2100e3;                % Altitude of perigee (m) 

inc = 28.5*(pi/180);        % Inclination (radians)
w = 0;                      % Argument of perigee (radians)
bigOmega = 0;               % Right ascension of the ascending node - RAAN (radians)
m0 = 0;                     % Mean anomaly at epoch (radians)

[t, pos, vel, nu, e] = getOrbit(hA, hP, inc, w, bigOmega, m0, dt);

t = 0:dt:nOrbits*t(end) + (nOrbits-1)*dt;
pos = repmat(pos, nOrbits, 1);
vel = repmat(vel, nOrbits, 1);
e = repmat(e, length(t), 1);
nu = repmat(nu, nOrbits, 1);
nudot = getNuDot(nu, dt);
sun = getSunEphemeris(date, dt, t(end));

% Initial quaternion
bz =  vel(1,:)/norm(vel(1,:));
by = cross(pos(1,:), vel(1,:))/norm(cross(pos(1,:),vel(1,:)));
bx = -cross(bz, by);
q0 = getQuaternion([bx; by; bz]);

% Initial angular velocities in body frame
w0 = [0.01 0.01 0.01];
torqueFuncs = {@getGravityGradient {e nu nudot},
               @getSolarRadiationPressure {pos sun}};
[q, w, tau] = propagateAttitude(sc, t, dt, q0, w0, torqueFuncs);

Lt = length(t);
d = zeros(Lt, 1);
d(1) = datenum(date);
for i = 2:Lt
    d(i) = addtodate(d(i-1), dt, 'second');
end
d = datevec(d);

% Ground stations
groundStations = getGroundStations();
gsFields = fieldnames(groundStations);
labels = {};
barAccessTimes = [];
for i = 1:numel(gsFields)
    gs = groundStations.(gsFields{i});
    lla = repmat(gs.LLA, length(d), 1);
    % Populate ground station position vector for date range
    gs.pos = lla2eci(lla, d);
    if gs.isSLR == 1
        gs = getSLRAccess(gs, retros, pos, vel, q);
        labels{end+1} = gs.name;
    end
    groundStations.(gsFields{i}) = gs;
end

% Rotations
hF1 = figure(1);
set(hF1, 'name', ['Rotations (' scName ')'], 'NumberTitle', 'off');
subplot(3,1,1);
plot(t/86400, w(:,1));
ylabel('\omega_1 (rads/sec)');
subplot(3,1,2);
plot(t/86400, w(:,2));
ylabel('\omega_2 (rads/sec)');
subplot(3,1,3);
plot(t/86400, w(:,3));
ylabel('\omega_3 (rads/sec)');
xlabel('Time (days)');

% Torques
[~, c] = size(tau);
n = c/3;

hF2 = figure(2);
set(hF2, 'name', ['External Torques (' scName ')'], 'NumberTitle', 'off');
subplot(3,1,1);
hold on;
for i = 1:n
    plot(t/86400, tau(:,(i*3-3)+1));
end
ylabel('M_1 (N m)');
subplot(3,1,2);
for i = 1:n
    plot(t/86400, tau(:,(i*3-3)+2));
    hold on;
end
ylabel('M_2 (N m)');
subplot(3,1,3);
for i = 1:n
    plot(t/86400, tau(:,(i*3-3)+2));
    hold on;
end
ylabel('M_3 (N m)');
xlabel('Time (days)');

% Ground station range and zenith angle
hF3 = figure(3);
set(hF3, 'name', 'Range and Zenith Angles', 'NumberTitle', 'off');
subplot(2,1,1);
hold on;
subplot(2,1,2);
hold on;
colors = {};
for i = 1:numel(gsFields)
    gs = groundStations.(gsFields{i});
    if isfield(gs, 'range')
        subplot(2,1,1);
        hRange = plot(t/86400, gs.range/1000);
        % Save color for future plotting
        groundStations.(gsFields{i}).color = hRange.Color;
        colors{end+1} = hRange.Color;
    end
    if isfield(gs, 'angle')
        subplot(2,1,2);
        plot(t/86400, gs.angle*(180/pi));
    end
end
subplot(2,1,1);
legend(labels);
ylabel('Slant Range (km)');
grid on;
subplot(2,1,2);
xlabel('Time (days)');
ylabel('Zenith Angle (degrees)');
grid on;

% Velocity aberrations
hF4 = figure(4);
set(hF4, 'name', ['Velocity Aberrations (' scName ')'], 'NumberTitle', 'off');
hold on;
for i = 1:numel(gsFields)
    gs = groundStations.(gsFields{i});
    if isfield(gs, 'alpha')
        scatter(t(gs.indices)/86400, gs.alpha*1e6, 'filled');
    end
end
grid on;
scatterLegend(colors, labels);
xlabel('Time (days)');
ylabel('Velocity Aberration (\murad)');

% Incidence angles
hF5 = figure(5);
set(hF5, 'name', ['Incidence Angles (' scRetrosName ')'], 'NumberTitle', 'off');
hold on;
for i = 1:numel(gsFields)
    gs = groundStations.(gsFields{i});
    if isfield(gs, 'thetaInc')
        c = makeColorMap([1 1 1], gs.color, 100);
        c = c(10:100,:);
        for j = 1:numel(retros.retros)
            retro = retros.retros(j);
            % Generate color gradient for points based on the effective area
            eta = getAreaFactor(retro.n, gs.thetaInc(:,j));
            bins = discretize(eta*90, 1:1:90);
            bins(isnan(bins)) = 1;
            csorted = c(bins,:);
            scatter(t(gs.indices)/86400, gs.thetaInc(:,j)*(180/pi), 25, csorted, 'filled');
        end
    end
end
grid on;
scatterLegend(colors, labels);
xlabel('Time (days)');
ylabel('Incidence Angle (degrees)');
colormap(c);
hcb = colorbar('Direction', 'reverse');
set(get(hcb, 'Title'), 'String', 'Effective Area');

% Optical cross-section
hF6 = figure(6);
set(hF6, 'name', ['Optical Cross-sections (' scRetrosName ')'], 'NumberTitle', 'off');
hold on;
for i = 1:numel(gsFields)
    gs = groundStations.(gsFields{i});
    if isfield(gs, 'sigma')
        scatter(t(gs.indices)/86400, gs.sigma, 'filled');
    end
end
set(gca, 'yscale', 'log');
grid on;
scatterLegend(colors, labels);
xlabel('Time (days)');
ylabel('Optical Cross-section (m^2)');

% Mean photoelectrons detected
hF7 = figure(7);
set(hF7, 'name', ['Mean Photoelectrons Detected (' scRetrosName ')'], 'NumberTitle', 'off');
hold on;
for i = 1:numel(gsFields)
    gs = groundStations.(gsFields{i});
    if isfield(gs, 'n')
        gs.n(gs.n==0) = NaN;
        % Generate color gradient for points based on the probability
        c = makeColorMap([1 1 1], gs.color, 100);
        c = c(10:100,:);
        bins = discretize(gs.Pdb*90, 1:1:90);
        bins(isnan(bins)) = 1;
        csorted = c(bins,:);
        scatter(t(gs.indices)/86400, gs.n, 25, csorted, 'filled');
    end
end
set(gca, 'yscale', 'log');
grid on;
scatterLegend(colors, labels);
xlabel('Time (days)');
ylabel('Mean Photoelectrons Detetected per Pulse');
colormap(c);
hcb = colorbar();
set(get(hcb, 'Title'), 'String', 'Probability of Detection');

% Probability of detection
hF8 = figure(8);
set(hF8, 'name', ['Probability of Detection (' scRetrosName ')'], 'NumberTitle', 'off');
hold on;
boxProb = [];
boxG = [];
for i = 1:numel(gsFields)
    gs = groundStations.(gsFields{i});
    if isfield(gs, 'Pdb')
        scatter(t(gs.indices)/86400, gs.Pdb, 'filled');
    end
end
ylim([0 1]);
grid on;
scatterLegend(colors, labels);
xlabel('Time (days)');
ylabel('Probability of Detection');

% Save all figures
if ~exist(outFolder, 'dir')
   mkdir(outFolder);
end
hF = [hF1 hF2 hF3 hF5 hF6 hF7 hF8];
for i = 1:length(hF)
    name = hF(i).Name;
    % Strip out s/c and retro config names
    k = strfind(name, ['(' scName]);
    if ~isempty(k)
        name = name(1:k-2);
    end
    % Strip out spaces
    filename = [outFolder '/' name(name~=' ') '.fig'];
    savefig(hF(i), filename, 'compact');
end
