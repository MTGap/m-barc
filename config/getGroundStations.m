function [groundStations] = getGroundStations()
    % The positions and properties of ground stations
    %
    % Michael Gapczynski
    %-------------------------------------------------
    
    % Ann Arbor, Michigan - SRB
    % 42.2943665 N 83.7118976 W 172 m
    gsAA = struct();
    gsAA.name = 'Ann Arbor';
    gsAA.LLA = [42.2943665 -83.7118976 172];
    gsAA.isSLR = 0;

    % Default SLR ground station properties
    gsSLR = struct();
    gsSLR.isSLR = 1;
    gsSLR.maxZenith = 70*(pi/180);    % Maximum Zenith Angle (radians)
    gsSLR.lambda = 532.1e-9;          % Laser wavelength (m)
    gsSLR.Epulse = 10e-3;             % Energy per pulse (J)
    gsSLR.thetaD = 20*4.84814e-6;     % Beam divergence half-angle (radians)
    gsSLR.deltaThetaP = 0;            % Beam pointing error (radians)
    gsSLR.deltaThetaJ = 2*4.84814e-6; % RMS tracking mount jitter (radians)
    gsSLR.nt = 0.85;                  % Transmitter optical throughput efficiency
    gsSLR.nr = 0.55;                  % Receiver optical throughput efficiency
    gsSLR.nc = 0.18;                  % Detector quantum efficiency
    gsSLR.Ar = 0.4;                   % Effective telescope receive area (m^2)
    
    % Zimmerwald, Switzerland (SLR)
    % http://ilrs.gsfc.nasa.gov/network/stations/active/ZIML_general.html
    % 46.8772 N 7.4652 E 951.2 m
    gsZ = gsSLR;
    gsZ.name = 'Zimmerwald';
    gsZ.LLA = [46.8772 7.4652 951.2];
    
    % Greenbelt, Maryland (SLR)
    % http://ilrs.gsfc.nasa.gov/network/stations/active/GODL_general.html
    % 39.0206 N 76.82770 W 19.184 m
    gsGr = gsSLR;
    gsGr.name = 'Greenbelt';
    gsGr.LLA = [39.0206 -76.82770 19.184];
    
    % McDonald, Texas (SLR)
    % http://ilrs.gsfc.nasa.gov/network/stations/active/MDOL_general.html
    % 30.6802 N 255.9848 E 2006.2210 m
    gsMc = gsSLR;
    gsMc.name = 'McDonald';
    gsMc.LLA = [30.6802 255.9848 2006.221];

    groundStations = struct('gsAA', gsAA, 'gsZ', gsZ, 'gsGr', gsGr, 'gsMc', gsMc);
end
