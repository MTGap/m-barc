close all;

addpath(genpath('config'));
addpath(genpath('lib'));

% Tethers Unlimited HYDROS Mission CONOPS
%
% Michael Gapczynski
%---------------------------------------- 

outFolder = 'output/hydros';

if ~exist(outFolder, 'dir')
   mkdir(outFolder);
end

% Record the animation
% vW = VideoWriter([outFolder '/orbit.avi']);
% vW.Quality = 100;
% open(vW);
vW = 0;

Re = 6378e3;                    % Radius of Earth (m)
AU = 1.4960e11;                 % Astronomical unit (m)

% Time
dt = 0.5;                       % Time step (s)
date = datetime(2018, 8, 20);   % Year, month, day
dateT = date;

% Solar array
solarArray = getSolarArrayZ();

% XB1
PXB1N = 6.3;                    % XB1 Nominal Power (W)
thetaDotXB1 = 10*(pi/180);      % Maneuver rate (rads/s)

% Battery
CBattery0 = 25;                 % Initial battery capacity (W hr)
CBatteryMax = CBattery0;        % Maximum battery capacity (W hr)
DOD = 0.6;                      % Depth of discharge

% HYDROS
PEz0 = 20;                      % HYDROS Electrolyzer power (W)
PEz = [];
VEz = 1.8;                      % Electrolyzer voltage (V)
MP = 0.01;                      % Mass of propellant pulse (kg)

tH2O = 132;                     % Electrolysis time (s)
tIH2O = ceil(tH2O/dt);

% Orbit Parameters
inc = 28.5*(pi/180);            % Inclination (radians)
w = 0;                          % Argument of perigee (radians)
bigOmega = 0;                   % Right ascension of the ascending node - RAAN (radians)
m0 = 0;                         % Mean anomaly at epoch (radians)

% Read altitude data from STK
data = xlsread('data/altitudes.xlsx');
hAV = 1e3*data(1:10, 4);
hPV = 1e3*data(1:10, 5);
% Duplicate rows to match up corresponding apogees and perigees
hAV = kron(hAV, ones(2,1));
hPV = [hPV(1); kron(hPV(2:end), ones(2,1))];

% Burns
nBA = 75;                       % Number of burns per apogee
nMA = 30;                       % Number of total maneuvers at apogee
tBA = 0.5;                      % Duration of burn at apogee (s)
nBP = 5;                        % Number of burns per perigee
nMP = 30;                       % Number of total maneuvers at perigee
tBP = 0.5;                      % Duration of burn per perigee (s)

tMA = nBA*tH2O + nBA*tBA;       % Duration of maneuver at apogee (s)
tMP = nBP*tH2O + nBP*tBP;       % Duration of maneuver at perigee (s)

tIMA = ceil(tMA/dt);
tIMP = ceil(tMP/dt);

t = [];
pos = [];
vel = [];
q = [];

% Calculate orbits from altitude data
for i = 1:length(hPV)
    if mod(i, 2) == 0
        % Start at apogee
        m0 = pi;
    else
        % Start at perigee
        m0 = 0;
    end
    [tT, posT, velT] = getOrbit(hAV(i), hPV(i), inc, w, bigOmega, m0, dt);
    tOrbit = length(tT);
    tHalf = ceil(tOrbit/2);
    % Concatenate orbits
    if i == 1
        t = tT(1:tHalf);
        bz =  velT(1,:)/norm(velT(1,:));
        by = cross(posT(1,:), velT(1,:))/norm(cross(posT(1,:),velT(1,:)));
        bx = -cross(bz, by);
        q(1,:) = getQuaternion([bx; by; bz]);
    else
        t = [t t(end) + dt + tT(1:tHalf)];    
    end
    pos = [pos; posT(1:tHalf,:)];
    vel = [vel; velT(1:tHalf,:)];
    % Determine attitude
    sun = getSunEphemeris(dateT, dt, tT(tHalf));
    dateT = dateT+seconds(tT(tHalf));
    % Burn counters
    nBAC = 0;
    nBPC = 0;
    for j = 1:length(sun)
        pointToSun = 0;
        % Maneuvers
        if mod(i, 2) == 0 && j < tIMA && nBAC <= nBA
            % Burn at apogee
            if dt <= tBA
                % Model individual burns if timestep allows
                if j >= nBAC*tIH2O
                    scV = [0 0 1];
                    pVI = -velT(j,:)/norm(velT(j,:));
                    PEz(end+1) = 0;
                    nBAC = nBAC+1;
                elseif j >= nBAC*tIH2O - 10/dt
                    % Orient s/c for upcoming burn
                    scV = [0 0 1];
                    pVI = -velT(j,:)/norm(velT(j,:));
                    PEz(end+1) = PEz0;
                else
                    PEz(end+1) = PEz0;
                    pointToSun = 1;
                end
            else
                scV = [0 0 1];
                pVI = -velT(j,:)/norm(velT(j,:));
                PEz(end+1) = PEz0;
            end
        elseif i > 1 && mod(i, 2) == 1 && j < tIMP && nBPC <= nBP
            % Burn at perigee
             if dt <= tBP
                % Model individual burns if timestep allows
                if j >= nBPC*tIH2O
                    scV = [0 0 1];
                    pVI = -velT(j,:)/norm(velT(j,:));
                    PEz(end+1) = 0;
                    nBPC = nBPC+1;
                elseif j >= nBPC*tIH2O - 10/dt
                    % Orient s/c for upcoming burn
                    scV = [0 0 1];
                    pVI = -velT(j,:)/norm(velT(j,:));
                    PEz(end+1) = PEz0;
                else
                    PEz(end+1) = PEz0;
                    pointToSun = 1;
                end
            else
                scV = [0 0 1];
                pVI = -velT(j,:)/norm(velT(j,:));
                PEz(end+1) = PEz0;
            end
        else
            pointToSun = 1;
            PEz(end+1) = 0;
        end
        if pointToSun == 1
            % Point solar array towards Sun when possible
            if ~inSunlight(pos(j,:), AU*sun(j,:))
                q(end+1,1:4) = q(end,1:4);
                continue;
            else
                scV = solarArray{1}.normal;
                pVI = -sun(j,:);
            end
        end
        pV = i2b(q(end,:), pVI);    % Vector to point s/c vector to
        angVel = thetaDotXB1*cross(pV, scV);
        angVelNorm = norm(angVel);
        co = cos(0.5*angVelNorm);
        si = sin(0.5*angVelNorm);
        n1 = angVel(1)/angVelNorm;
        n2 = angVel(2)/angVelNorm;
        n3 = angVel(3)/angVelNorm;
        qw1 = n1*si;
        qw2 = n2*si;
        qw3 = n3*si;
        qw4 = co;
        om = [qw4  qw3 -qw2 qw1;
           -qw3  qw4  qw1 qw2;
            qw2 -qw1  qw4 qw3;
           -qw1 -qw2 -qw3 qw4];
        q(end+1,1:4) = (om*q(end,1:4)')';
    end
end

q = q(1:end-1,1:4);

[PSolar, solarArray] = getSolarPower(solarArray, date, dt, t(end), pos, q);

% Solar panel incidence angles
hF1 = figure(1);
set(hF1, 'name', 'Solar Panel Incidence Angles', 'NumberTitle', 'off');
hold on;
labels = [];
for i = 1:numel(solarArray)
    sp = solarArray{i};
    plot(t/86400, sp.thetaInc*(180/pi));
    labels{end+1} = sp.name;
end
legend(labels);
xlabel('Time (days)');
ylabel('Sun Incidence Angle (degrees)');

% Power consumption
hF2 = figure(2);
set(hF2, 'name', 'Power Consumption', 'NumberTitle', 'off');
PXB1 = PXB1N*ones(1, length(t));
plot(t/86400, -PXB1);
hold on;
plot(t/86400, -PEz);
plot(t/86400, PSolar);
ylim([-30 30]);
legend('XB1', 'Electrolyzer', 'Solar Array');
xlabel('Time (days)');
ylabel('Power (W)');

% Battery capacity
hF3 = figure(3);
set(hF3, 'name', 'Battery Capacity', 'NumberTitle', 'off');
PConsumed = PXB1+PEz;
CBattery = zeros(length(t), 1);
CBattery(1) = CBattery0;

% Battery charging state
% 0 = Discharging, 1 = Charging
isCharging = 0;
cycles = 0;
for i = 2:length(t)
   PDiff = (PSolar(i) - PConsumed(i))*(dt/3600);
   if PDiff < 0
       % Discharge
       CBattery(i) = CBattery(i-1) + PDiff;
   elseif isCharging
       % Charge
       if CBattery(i-1) + PDiff < CBatteryMax
           CBattery(i) = CBattery(i-1) + PDiff;
       else
           CBattery(i) = CBatteryMax;
       end
   else
       CBattery(i) = CBattery(i-1);
   end
   if CBattery(i) <= CBatteryMax*(1-DOD)
       % Start charging when depth of discharge is reached
       isCharging = 1;
   elseif CBattery(i) == CBatteryMax && isCharging
       % Stop charging when max capacity is reached
       isCharging = 0;
       cycles = cycles+1;
   end
end
plot(t/86400, CBattery);
xlabel('Time (days)');
ylabel('Battery Capacity (W hr)');
ylim([0 30]);

% Save all figures
hF = [hF1 hF2 hF3];
for i = 1:length(hF)
    name = hF(i).Name;
    % Strip out spaces
    filename = [outFolder '/' name(name~=' ') '.fig'];
    savefig(hF(i), filename, 'compact');
end

% Plot the orbit
hO = figure(4);
posEnd = length(pos)-tOrbit;
plot3(pos(1:posEnd,1), pos(1:posEnd,2), pos(1:posEnd,3));
hold on;
% Plot the last orbit as a different color
plot3(pos(posEnd:length(pos),1), pos(posEnd:length(pos),2), pos(posEnd:length(pos),3));
axis equal;
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
rotEdot = diff(rotE)*dt;
tEnd = length(pos);
rotE = mod(rotE(1):rotEdot:rotE(1)+rotEdot*tEnd, 360);
rotate(hMesh, [0 0 1], rotE(1), [0 0 0]);

x = pos(1,1);
y = pos(1,2);
z = pos(1,3);

% s/c
hSC = plot3(x, y, z, 'ob', 'MarkerSize', 6, 'MarkerFaceColor', 'b');

% Body frame
bA = getAttitude(q(1,:));
ax = [x y z; 7e6*bA(1,:) + [x y z]];
ay = [x y z; 7e6*bA(2,:) + [x y z]];
az = [x y z; 7e6*bA(3,:) + [x y z]];
hBx = plot3(ax(:,1), ax(:,2), ax(:,3), 'r');
hBy = plot3(ay(:,1), ay(:,2), ay(:,3), 'g');
hBz = plot3(az(:,1), az(:,2), az(:,3), 'b');

% Sun vector
sun = getSunEphemeris(date, dt, t(end));
as = [0 0 0; 1e7*sun(1,:)];
hSun = plot3(as(:,1), as(:,2), as(:,3), 'Color', [0.85 0.325 0.098], 'LineWidth', 2);

% ECEF frame
ax = [0 0 0; 1e7*cosd(rotE(1)) 1e7*sind(rotE(1)) 0];
ay = [0 0 0; -1e7*sind(rotE(1)) 1e7*cosd(rotE(1)) 0];
az = [0 0 0; 0 0 1e9];
hEX = plot3(ax(:,1), ax(:,2), ax(:,3), 'k', 'LineWidth', 2);
hEY = plot3(ay(:,1), ay(:,2), ay(:,3), 'k', 'LineWidth', 2);
hEZ = plot3(az(:,1), az(:,2), az(:,3), 'k', 'LineWidth', 2);

handles = struct('hSC', hSC, ...
                 'hBx', hBx, 'hBy', hBy, 'hBz', hBz, ...
                 'hSun', hSun, 'hEX', hEX, 'hEY', hEY, 'hMesh', hMesh);

% Animate the s/c position and attitude
at = timer('ExecutionMode', 'fixedRate', 'Period', 0.05, ...
          'TimerFcn', {@animateSC, vW, tEnd, handles, pos, q, sun, rotE});
start(at);

% Stop animation timer if the figure is closed
hO.DeleteFcn = @(~, ~) stop(at);
