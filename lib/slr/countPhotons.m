function [n] = countPhotons(gs, range, sigma, Ta, Tc)
    % Count the number of photons returned by a retroreflector
    % 
    % gs: Ground station structure array
    %   - gs.lambda: Laser wavelength (m)
    %   - gs.Epulse: Energy per pulse (J)
    %   - gs.thetaD: Beam divergence half-angle (radians)
    %   - gs.deltaThetaP: Beam pointing error (radians)
    %   - gs.deltaThetaJ: RMS tracking mount jitter (radians)
    %   - gs.nt: Transmitter optical throughput efficiency
    %   - gs.nr: Receiver optical throughput efficiency
    %   - gs.nc: Detector quantum efficiency
    %   - gs.Ar: Effective telescope receive area (m^2)
    % sigma: Optical cross-section (m^2)
    % Ta: One-way atmospheric transmission
    % Tc: One-way cirrus transmission
    %
    % countPhotons(gs.Z, 2100e3, 1e6, 0.8, 1)
    %
    % Michael Gapczynski
    %---------------------------------------------------------
    
    h = 6.62607e-34;        % Planck's constant (J s)
    c = 299792458;          % Speed of light (m/s)

    v = c/gs.lambda;        % Frequency (Hz)

    % Transmitter gain
    Gt = (8/gs.thetaD^2)*exp(-2*(gs.deltaThetaP/gs.thetaD)^2);

    % Range loss
    Ls = (1./(4*pi*range.^2)).^2;
    
    % Tracking loss
    Lt = 1/(1+(gs.deltaThetaJ/gs.thetaD)^2);

    n = (gs.Epulse./(h.*v)).*gs.nt.*Gt.*Lt.*sigma.*Ls*gs.Ar*gs.nr*gs.nc.*Ta.^2.*Tc.^2;
end
