close all;

addpath(genpath('lib'));

% Test conversion between text and ASCII bitmap
%
% Michael Gapczynski
%----------------------------------------------

textFile = 'data/sample.txt';
outFolder = 'output/payload';

% 100 nm x 100 nm etched squares
squareSL = 100e-9;    % Side length of etched square (m)
chipSL = 1*0.0254;    % Side length of payload chip (m)

cols = chipSL/squareSL;
maxCols = cols-mod(cols, 7);

text = fileread(textFile);
bitmap = convertTextToASCIIBitmap(text, maxCols);
imshow(bitmap);

textOut = convertASCIIBitmapToText(bitmap);
disp(textOut);

if ~exist(outFolder, 'dir')
   mkdir(outFolder);
end

% Determine original filename
[~, name, ~] = fileparts(textFile);

% Save bitmap
imwrite(bitmap, [outFolder '/' name '.bmp']);

% Save converted text
fid = fopen([outFolder '/' name 'Out.txt'], 'wt');
fprintf(fid, '%s\n', textOut);
fclose(fid);

