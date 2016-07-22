function [solarArray] = getSolarArrayZ()
    % A 4 panel 2.5U solar array configuration on the -Z face
    %
    % Michael Gapczynski
    %---------------------------------------------------------
    
    area = 0.0634;     % Area (m^2)
    eta = 0.3;         % Efficiency
    
    sp = struct('area', area, 'eta', eta);
    
    % -Z
    sp1 = sp;
    sp1.name = '-Z';
    sp1.normal = [0 0 -1];

    solarArray = {sp1};
end
