function [] = scatterLegend(colors, labels)
    % Create a legend for scatter plots
    %
    % colors: 1xn cell array
    % labels: 1xn cell array
    %
    % MATLAB's legend function oftentimes displays incorrectly
    % for scatter plots. This function creates a fake, hidden 
    % plot to attach the legend to instead. This seems to correct
    % most legend problems.
    %
    % Michael Gapczynski
    %-------------------------------------------------------------

    if length(colors) ~= length(labels)
        error('The length of colors and labels must be the same');
    end
    n = length(colors);
    hFake = plot(NaN(2,n));
    for i = 1:n
        hFake(i).Marker = 'o';
        hFake(i).MarkerFaceColor = colors{i};
        hFake(i).MarkerEdgeColor = colors{i};
        hFake(i).LineStyle = 'none';
    end
    legend(hFake, labels);
end
