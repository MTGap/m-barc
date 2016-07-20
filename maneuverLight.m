close all;

addpath(genpath('lib'));

% Plot the parameter space of deployment date and RAAN
% for eclipse times during maneuvers
%
% Michael Gapczynski
%-----------------------------------------------------

outFolder = 'output/light';

AU = 1.4960e11;         % Astronomical unit (m)

% Date range of plots
dateRange = datenum({'01-Jan-2018 00:00:00'; '31-Dec-2018 23:00:00'});
dates = datevec(dateRange(1):1:dateRange(2));

% Orbit Parameters
hA = 35786e3;           % Altitude of apogee (m)
hP = 200e3;             % Altitude of perigee (m)

inc = 28.5*(pi/180);    % Inclination (radians)
w = 180*(pi/180);       % Argument of perigee (typically 0 or 180 for GTO) (radians)
bigOmega = 0:0.01:2*pi; % Right ascension of the ascending node - RAAN (radians)
m0 = 0;                 % Mean anomaly at epoch (radians)

dt = 10;                % Orbit time step (s)

% HYDROS
tH2O = 132;     % Electrolysis time (s)

% Burns
nBA = 75;       % Number of burns per apogee
nMA = 30;       % Number of total maneuvers at apogee
tBA = 0.5;      % Duration of burn at apogee (s)
nBP = 5;        % Number of burns per perigee
nMP = 30;       % Number of total maneuvers at perigee
tBP = 0.5;      % Duration of burn per perigee (s)

tMA = nBA*tH2O + nBA*tBA;   % Duration of maneuver at apogee (s)
tMP = nBP*tH2O + nBP*tBP;   % Duration of maneuver at perigee (s)

% Orbit time step index lengths
tIMA = ceil(tMA/dt);        
tIMP = ceil(tMP/dt);

Ldates = length(dates);
LbigOmega = length(bigOmega);

pos = cell(Ldates, 1);
tHalf = zeros(Ldates, 1);

% Compute orbits
for i = 1:length(bigOmega)
    [tT, posT, ~] = getOrbit(hA, hP, inc, w, bigOmega(i), m0, dt);
    pos{i} = posT; 
    tHalf(i) = ceil(length(tT)/2);
end

lightA = zeros(LbigOmega, Ldates);
lightP = zeros(LbigOmega, Ldates);

for i = 1:Ldates
    date = datetime(dates(i,:));
    sun = getSunEphemeris(date, 1, 1);
    sun = sun(1,:);
    for j = 1:LbigOmega
        posT = pos{j};
        tHalfT = tHalf(j);
        % Calculate average sunlight duration for maneuvers
        lightA(j, i) = mean(inSunlight(posT(tHalfT:tHalfT+tIMA,:), AU*repmat(sun, tIMA+1, 1)));
        lightP(j, i) = mean(inSunlight(posT(1:1+tIMP,:), AU*repmat(sun, tIMP+1, 1)));
    end   
end

% Apogee Maneuver Eclipse
hF1 = figure(1);
set(hF1, 'name', 'Apogee Maneuver Eclipse', 'NumberTitle', 'off');
contourf(1:Ldates, bigOmega*(180/pi), lightA, 100, 'EdgeColor', 'none');
grid on;
colormap bone;
ylabel('RAAN (degrees)');
xlabel('Deployment Date');
xlim([0 length(dates)]);
ylim([bigOmega(1)*(180/pi) bigOmega(end)*(180/pi)]);
datetick('x', 'keeplimits');

% Perigee Maneuver Eclipse
hF2 = figure(2);
set(hF2, 'name', 'Perigee Maneuver Eclipse', 'NumberTitle', 'off');
contourf(1:Ldates, bigOmega*(180/pi), lightP, 100, 'EdgeColor', 'none');
grid on;
colormap bone;
xlabel('Deployment Date');
ylabel('RAAN (degrees)');
xlim([0 length(dates)]);
ylim([bigOmega(1)*(180/pi) bigOmega(end)*(180/pi)]);
datetick('x', 'keeplimits');

% Apogee and Perigee Maneuver Eclipse
hF3 = figure(3);
set(hF3, 'name', 'Maneuver Eclipse', 'NumberTitle', 'off');
% Scale colors based on number of burns
lights = (nBA/(nBA+nBP))*lightA+(nBP/(nBA+nBP))*lightP;
contourf(1:Ldates, bigOmega*(180/pi), lights, 100, 'EdgeColor', 'none');
grid on;
colormap bone;
xlabel('Deployment Date');
ylabel('RAAN (degrees)');
xlim([0 length(dates)]);
ylim([bigOmega(1)*(180/pi) bigOmega(end)*(180/pi)]);
datetick('x', 'keeplimits');

% Save all figures
if ~exist(outFolder, 'dir')
   mkdir(outFolder);
end
hF = [hF1 hF2 hF3];
for i = 1:length(hF)
    name = hF(i).Name;
    % Strip out spaces
    filename = [outFolder '/' name(name~=' ') '.fig'];
    savefig(hF(i), filename, 'compact');
end

