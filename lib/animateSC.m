function animateSC(timer, ~, vW, tEnd, handles, pos, q, sun, rotE)
    % Animate the position and attitude of the s/c and other objects
    %
    % timer: Animation timer (timer)
    % vW: VideoWriter
    % tEnd: Index of end of animation
    % handles: Figure handles to the objects to animate
    %   - hSC: s/c
    %   - hBx, hBy, hBz: s/c body frame vectors
    %   - hSun: Sun vector
    %   - hEX, hEY: ECEF frame vectors
    %   - hMesh: Earth map mesh
    %   - groundStations: Ground station struct
    %
    % pos: s/c position vector in ECI frame
    % q: s/c attitude quaternion
    % sun: Sun unit vector in ECI frame
    % rotE: Rotation of Earth per time step (degrees)
    %
    % Michael Gapczynski
    %----------------------------------------------------------------------
    
    AU = 1.4960e11;         % Astronomical unit (m)
    
    i = mod(timer.TasksExecuted, tEnd) + 1;
    
    x = pos(i,1);
    y = pos(i,2);
    z = pos(i,3);
    
    % s/c
    set(handles.hSC, 'XData', x, 'YData', y, 'ZData', z); 
    if inSunlight(pos(i,:), AU*sun(i,:))
        set(handles.hSC, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    else
        set(handles.hSC, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    end

    % Body frame
    bA = getAttitude(q(i,:));
    ax = [x y z; 7e6*bA(1,:) + [x y z]];
    ay = [x y z; 7e6*bA(2,:) + [x y z]];
    az = [x y z; 7e6*bA(3,:) + [x y z]];
    set(handles.hBx, 'XData', ax(:,1), 'YData', ax(:,2), 'ZData', ax(:,3));
    set(handles.hBy, 'XData', ay(:,1), 'YData', ay(:,2), 'ZData', ay(:,3));
    set(handles.hBz, 'XData', az(:,1), 'YData', az(:,2), 'ZData', az(:,3));

    % Sun vector
    as = [0 0 0; 1e7*sun(i,:)];
    set(handles.hSun, 'XData', as(:,1), 'YData', as(:,2), 'ZData', as(:,3));
    
    % ECEF frame
    ax = [0 0 0; 1e7*cosd(rotE(i)) 1e7*sind(rotE(i)) 0];
    ay = [0 0 0; -1e7*sind(rotE(i)) 1e7*cosd(rotE(i)) 0];
    set(handles.hEX, 'XData', ax(:,1), 'YData', ax(:,2), 'ZData', ax(:,3));
    set(handles.hEY, 'XData', ay(:,1), 'YData', ay(:,2), 'ZData', ay(:,3));
    
    % Earth
    if i == 1 && timer.TasksExecuted > 2
        % Reset Earth to original rotation
        alpha = rotE(tEnd) - rotE(1);
        rotate(handles.hMesh, [0 0 1], -alpha, [0 0 0]);
    else
        rotEdot = diff(rotE(1:2));
        rotate(handles.hMesh, [0 0 1], rotEdot, [0 0 0]);
    end
    
    % Ground stations
    if isfield(handles, 'groundStations')
        gsFields = fieldnames(handles.groundStations);
        for j = 1:numel(gsFields)
           gs = handles.groundStations.(gsFields{j});
           posGS = gs.pos(i,:);
           set(gs.hGS, 'XData', posGS(1), 'YData', posGS(2), 'ZData', posGS(3));
           if gs.isSLR == 1
               % Turn on laser if in accessible region
               if gs.range(i) <= 10000e3 && gs.angle(i) <= 1.2217
                   aLaser = [gs.pos(i,:); [x y z]];
                   set(gs.hLaser, 'LineStyle', '-');
                   set(gs.hLaser, 'XData', aLaser(:,1), 'YData', aLaser(:,2), 'ZData', aLaser(:,3));
               else
                   set(gs.hLaser, 'LineStyle', 'none');
               end
           end
        end
    end
    
    % Save frame to video
    if isa(vW, 'VideoWriter')
        if timer.TasksExecuted <= tEnd
            writeVideo(vW, getframe(gca));
        elseif timer.TasksExecuted == tEnd+1
            close(vW);
        end
    end
end
