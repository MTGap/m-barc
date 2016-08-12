close all;

addpath(genpath('lib'));

% Earth's magnetic field strength calculations for
% an orbit using the IGRF-12 magnetic field model
%
% http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
%
% Michael Gapczynski
%----------------------------------------------------

outFolder = 'output/magnetic';

Re = 6378e3;                    % Radius of Earth (m)

% Time
dt = 100;                       % Time step (s)
date = datetime(2018, 1, 20);   % Year, month, day

% Magnetorquer
D = 0.12;                       % Max dipole moment (A m^2)

% Orbit Parameters
hA = 35786e3;                   % Altitude of apogee (m)
hP = 200e3;                     % Altitude of perigee (m)

inc = 0;                        % Inclination (radians)
w = 0;                          % Argument of perigee (radians)
bigOmega = 0;                   % Right ascension of the ascending node - RAAN (radians)
m0 = pi;                        % Mean anomaly at epoch (radians)

% When the inclination is set to 0 degrees, the magnetic field strength
% in the North direction of the NED frame will match the magnetic field
% strength in the Z direction of the ECI frame.

[t, pos, vel, nu] = getOrbit(hA, hP, inc, w, bigOmega, m0, dt);

% Initialize date vector
Lt = length(t);
d = zeros(Lt, 1);
d(1) = datenum(date);
for i = 2:Lt
    d(i) = addtodate(d(i-1), dt, 'second');
end
d = datevec(d);

lat = zeros(Lt, 1);
lon = zeros(Lt, 1);
alt = zeros(Lt, 1);
BNED = zeros(Lt, 3);
BECI = zeros(Lt, 3);
for i = 1:length(d)
    lla = eci2lla(pos(i,:), d(i,:));
    lat(i) = lla(1);
    lon(i) = lla(2);
    alt(i) = lla(3);
    % Get components of the Earth's magnetic field
    B = igrf(datenum(d(i,:)), lat(i), lon(i), alt(i)/1000);
    % Convert from nT to T
    B = B/1e9;
    % Convert from NED to ECEF to ECI
    DCMNED = [-sind(lat(i))*cosd(lon(i)) -sind(lat(i))*sind(lon(i)) cosd(lat(i));
              -sind(lon(i)) cosd(lon(i)) 0;
              -cosd(lat(i))*cosd(lon(i)) -cosd(lat(i))*sind(lon(i)) -sind(lat(i))];
    DCMECI = dcmeci2ecef('IAU-2000/2006', d(i,:));
    BNED(i,1:3) = B;
    BECI(i,1:3) = DCMECI'*DCMNED'*B';
end

hF1 = figure(1);
set(hF1, 'name', 'Magnetic Field in NED Frame', 'NumberTitle', 'off');
plot(t/86400, BNED(:,1));
hold on;
plot(t/86400, BNED(:,2));
plot(t/86400, BNED(:,3));
legend('B_{North}', 'B_{East}', 'B_{Down}');
xlabel('Time (days)');
ylabel('Magnetic Field Strength (T)');

hF2 = figure(2);
set(hF2, 'name', 'Magnetic Field in ECI Frame', 'NumberTitle', 'off');
plot(t/86400, BECI(:,1));
hold on;
plot(t/86400, BECI(:,2));
plot(t/86400, BECI(:,3));
legend('B_x', 'B_y', 'B_z');
xlabel('Time (days)');
ylabel('Magnetic Field Strength (T)');

% Calculate ideal case magnetorquer torques
% with 3 magnetorquers on the NED axes
mag = zeros(Lt, 3);
magX = mag;
magX(:,1) = D;
magY = mag;
magY(:,2) = D;
magZ = mag;
magZ(:,3) = D;

tauX = cross(magX, BNED);
tauY = cross(magY, BNED);
tauZ = cross(magZ, BNED);

tau = tauX+tauY+tauZ;

hF3 = figure(3);
set(hF3, 'name', 'Magnetorquers Aligned with NED Frame', 'NumberTitle', 'off');
semilogy(t/86400, abs(tau(:,1)));
hold on;
semilogy(t/86400, abs(tau(:,2)));
semilogy(t/86400, abs(tau(:,3)));
legend('\tau_{North}', '\tau_{East}', '\tau_{Down}');
xlabel('Time (days)');
ylabel('Magnetorquer Torque (N m)');

hO = figure(4);
set(hO, 'name', 'Earth Magnetic Field', 'NumberTitle', 'off');
plot3(pos(:,1), pos(:,2), pos(:,3));
hold on;
axis(1.3*axis);
xlabel('x');
ylabel('y');
zlabel('z');

% Plot Earth
[XS, YS, ZS] = sphere(30);
hEarth = surf(XS*Re, YS*Re, ZS*Re, 'FaceColor', 'none', 'EdgeColor', 'none');
load topo;
axesm('globe', 'Geoid', Re);
hMesh = meshm(topo, topolegend);
demcmap(topo);

% Calculate rotation of Earth
JD = juliandate([date date+seconds(1)]);
rotE = JD2GMST(JD);
rotate(hMesh, [0 0 1], rotE(1), [0 0 0]);

latStart = 30:15:60;
lonStart = 0:30:330;
altStart = 0;
distance = -sign(latStart).*[30e3 70e3 150e3];
nsteps = abs(distance)/10;

% Calculate magnetic field lines
lat = zeros(max(nsteps(:))+1, numel(latStart)*numel(lonStart));
lon = zeros(max(nsteps(:))+1, numel(latStart)*numel(lonStart));
alt = zeros(max(nsteps(:))+1, numel(latStart)*numel(lonStart));
for index1 = 1:numel(latStart)
    for index2 = 1:numel(lonStart)
        [lat(1:nsteps(index1)+1, ...
            index1*(numel(lonStart)-1)+index2), lon(1:nsteps(index1)+1, ...
            index1*(numel(lonStart)-1)+index2), alt(1:nsteps(index1)+1, ...
            index1*(numel(lonStart)-1)+index2)] = ...
            igrfline(datenum(date), latStart(index1), lonStart(index2), ...
            altStart, 'geod', distance(index1), nsteps(index1));
    end
end
plot3m(lat, lon, alt*1000, 'r');
axis equal;

% Save all figures
if ~exist(outFolder, 'dir')
   mkdir(outFolder);
end
hF = [hF1 hF2 hF3 hO];
for i = 1:length(hF)
    name = hF(i).Name;
    % Strip out spaces
    filename = [outFolder '/' name(name~=' ') '.fig'];
    savefig(hF(i), filename, 'compact');
end
