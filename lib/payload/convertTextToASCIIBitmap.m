function [bitmap] = convertTextToASCIIBitmap(text, maxCols)
    % Convert text to an ASCII character encoded bitmap
    %
    % text: String of text to convert
    % maxCols: Maximum number of columns of the bitmap
    %          (has to be divisible by 7) 
    %
    % convertTextToASCIIBitmap(fileread('transcript.txt'), 49) 
    %
    % Michael Gapczynski
    %----------------------------------------------------------
    
    nBits = 7;
    
    if mod(maxCols, nBits) ~= 0
        error('Maximum number of columns has to be divisible by 7');
    end
    
    % Convert to binary
    dASCII = double(text);
    lookupTable = cell2mat(arrayfun(@(i)bitget((0:2^nBits - 1)', nBits+1-i), 1:nBits, 'UniformOutput', 0));
    bASCII = lookupTable(dASCII(:)+1,:);
 
    % Calculate the size to make the bitmap as square as possible
    n = numel(bASCII);
    sqn = ceil(sqrt(n));
    % Write left-to-right with whole 7-bit units in each row
    cols = sqn-mod(sqn, 7);
    % Don't exceed the specified maximum number of columns
    if cols > maxCols
       cols = maxCols; 
    end
    rows = ceil(n/cols);
    pad = ((rows*cols)-n)/7;
    % Add 0s to end in order to fill the last row of the bitmap
    bPASCII = padarray(bASCII, [pad 0], 'post');

    % Reshape into bitmap
    bitmap = reshape(bPASCII', cols, rows)';
    
    % Convert to logical to make it a binary image
    bitmap = logical(bitmap);
end

