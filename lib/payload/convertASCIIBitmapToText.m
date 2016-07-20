function [text] = convertASCIIBitmapToText(bitmap)
    % Convert an ASCII character encoded bitmap to text
    %
    % bitmap: ASCII character encoded bitmap
    %
    % Reads left-to-right from the top-left
    % Bitmap has to be formatted with whole 7-bit units in each row
    %
    % convertASCIIBitmapToText(bitmap) 
    %
    % Michael Gapczynski
    %---------------------------------------------------------------
    
    nBits = 7;
    
    [~, cols] = size(bitmap);
    if mod(cols, nBits) ~=0
        error('The bitmap is not properly formatted with whole 7-bit units in each row');
    end
    
    % Reshape into binary row vectors
    bASCII = reshape(bitmap', nBits, numel(bitmap)/nBits)';
    
    % Remove any extra 0s in the last row of the bitmap
    lastRow = find(sum(bASCII, 2), 1, 'last');
    bASCII = bASCII(1:lastRow, :);
    
    % Convert to characters
    dASCII = zeros(1, lastRow);
    base2 = 2.^(nBits-1:-1:0);
    for i = 1:lastRow
        dASCII(i) = sum(bASCII(i,:).*base2);
    end
    text = char(dASCII);
end

