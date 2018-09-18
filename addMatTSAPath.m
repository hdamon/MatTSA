% Configure Matlab Path to Include the full crlEEG
function addCRLEEGPath()
[currDir,~,~] = fileparts(mfilename('fullpath'));

addpath(currDir);
addpath([currDir '/external/wavelet']);

