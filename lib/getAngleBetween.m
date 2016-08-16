function [theta] = getAngleBetween(a, b)
    % Calculate the angle between two vectors
    %
    % a: Vector 1 (nx3)
    % b: Vector 2 (nx3)
    %
    % Michael Gapczynski
    %---------------------------------------------
    
    theta = atan2(vecnorm(cross(a, b, 2)), dot(a, b, 2));
end
