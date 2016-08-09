function [nudot] = getNuDot(nu, dt)
	% Calculate the rate of change of true anomaly
   	%
   	% nu: True anomaly (radians)
   	% dt: Time step (s)
   	%
    % Michael Gapczynski
    %-----------------------------------------------

    nudot = diff(nu)/dt;
    ids = 1:length(nudot);
    % Find and remove outliers that are resultant of vertical asymptotes
    outIds = abs(nudot - median(nudot)) > 3*std(nudot);
    nudot(outIds) = interp1(ids(~outIds), nudot(~outIds), ids(outIds), 'spline');
end
