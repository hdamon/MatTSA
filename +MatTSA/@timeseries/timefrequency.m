function out =  timefrequency(tseries,varargin)
%% Compute time-frequency decompositions of MatTSA.timeseries objects
%
% Inputs
% ------
%   tseries : A MatTSA.timeseries object to decompose
%   method  : Type of decomposition to compute
%               Valid values:
%                 'multitaper'  : Use pmtm for multitaper decomposition
%                 'spectrogram' : Use Matlab's spectrogram functionality.
%                 'fft'         : Use a fast fourier transform
%                 'eeglab'      : Use EEGlab's timefreq() function
%
% 'fft' and 'eeglab' are likely to be deprecated soon.
%
% Parameters for MultiTaper
% -------------------------
%    windowSize : Size of the time window (DEFAULT: 1024)
%            nw : Multitaper Parameter (DEFAULT: 3)
%     FFTLength : Length of FFT (DEFAULT: 1024)
%         freqs : Frequencies to compute decomposition at.
%       nOutput : Number of times to output at
%                   If nOutput<=1: Treated as a fraction of the total
%                                   samples
%                   Otherwise: Treated as an explicit number of samples
%
% MULTITAPER REQUIRES THE SIGNAL PROCESSING TOOLBOX
%
% Parameters for Spectrogram:
% ---------------------------
%   'window' :  Window to use for FFT computation 
%                 DEFAULT: hamming(2048)
%  'overlap' :  Number of samples to overlap the windows
%                 DEFAULT: 2038
%    'freqs' :  Vector defining the frequencies of the output
%
% Output
% ------
%    out.tfX : Time-frequency decomposition values
%    out.tx  : Times (in seconds) the decomposition was computed at
%    out.fx  : Frequencies the decomposition was computed at
%
% Part of the crlEEG project
% 2009-2018
%

%% Input Parsing
validTypes = {'multitaper' 'eeglab' 'spectrogram' 'wavelet'};
p = inputParser;
p.KeepUnmatched = true;
p.addOptional('method','wavelet',@(x) ismember(x,validTypes));
p.addRequired('tseries',@(x) isa(x,'MatTSA.timeseries'));
p.parse(tseries,varargin{:});

switch p.Results.method
  case 'fft'
    out = fft(tseries,p.Unmatched);
  case 'spectrogram'
    out = runspectrogram(tseries,p.Unmatched);
  case 'multitaper'
    out = multitaper(tseries,p.Unmatched);
  case 'eeglab'
    out = eeglab(tseries,p.Unmatched);
  case 'wavelet'
    out = continuousWaveletDecomp(tseries,p.Unmatched);
  otherwise
    error('Unknown decomposition type');
end;

end

%% Spectrogram Based Decomposition
function out = runspectrogram(tseries,varargin)
% Compute a time-frequency decomposition using Matlab's spectrogram
% 
% function out = runspectrogram(tseries,varargin)
%
% Inputs
% ------
%   tseries : MatTSA.timeseries object
%
% Param-Value Pairs
% --------
%   'window' :  Window to use for FFT computation 
%                 DEFAULT: hamming(2048)
%  'overlap' :  Number of samples to overlap the windows
%                 DEFAULT: 2038
%  'freqs'   :  Vector defining the frequencies of the output
%  
% REQUIRES THE SIGNAL PROCESSING TOOLBOX
% 

%% Input Parsing
p = inputParser;
p.addParameter('window',hamming(2048));
p.addParameter('overlap',2038); % every 10th sample
p.addParameter('freqs',linspace(0,50,100));
p.parse(varargin{:});

%% Computation
dataChans = logical(tseries.isChannelType('data'));
data = tseries.data(:,dataChans);
for i = 1:size(data,2)  
  [s(:,:,i),f,t] = spectrogram(data(:,i),...
    p.Results.window,...
    p.Results.overlap,...
    p.Results.freqs,...
    tseries.sampleRate);
end;

%% Output Parsing
out = MatTSA.tfDecomp('spectrogram',s,t+tseries.tRange(1),f,tseries.chanLabels(dataChans));
out.params = varargin;

end

function out = continuousWaveletDecomp(tseries,varargin)
% Compute a time-frequency decomposition using the continuous wavelet
% transform
%
%
% Inputs
% ------
%    tseries : A MatTSA.timeseries object
%
% Param-Value Pairs
% -----------------
%    'pad' :
%     'dj' :
%  '   s0' :
%     'j1' :
% 'mother' :
%  'param' :
%  'freqs' :
%
% This is computed using the function wavelet(), available from MatlabCentral
% at:
% https://www.mathworks.com/matlabcentral/fileexchange/20821


p = inputParser;
p.addRequired('tseries',@(x) isa(x,'MatTSA.timeseries'));
p.addParameter('pad',1);
p.addParameter('dj',0.1);
p.addParameter('s0',-1); 
p.addParameter('j1',-1);
p.addParameter('mother',-1);
p.addParameter('param',-1);
p.addParameter('freqs',[]);
p.parse(tseries,varargin{:});

PAD = p.Results.pad;
DJ = p.Results.dj;
S0 = p.Results.s0;
J1 = p.Results.j1;
MOTHER = p.Results.mother;
PARAM = p.Results.param;

DT = 1./tseries.sampleRate;

dataChans = logical(tseries.isChannelType('data'));
data = tseries.data(:,dataChans);

% Defaults from wavelet()
if (S0 ==-1 ), S0 = 2*DT; end;
if (DJ ==-1 ), DJ = 1/4; end;
nOut = fix((log(size(data,1)*DT/S0)/log(2))/DJ)+1;

doInterp = false;
if ~isempty(p.Results.freqs)        
  doInterp = true;
  nOut = numel(p.Results.freqs);
end

WAVE = zeros([nOut size(data)]);
for i = 1:size(data,2)  
%  [WAVE(:,:,i),PERIOD,SCALE,COI] = wavelet(data(:,i),DT,PAD,DJ,S0,J1,MOTHER,PARAM);
  [tmpWAVE,PERIOD,SCALE,COI] = wavelet(data(:,i),DT,PAD,DJ,S0,J1,MOTHER,PARAM);  
  
  tmpWAVE = flip(tmpWAVE,1);
  F = flip(1./PERIOD);
  Fout = F;  
  
  if doInterp
    % Interpolate, if needed
    Fout = p.Results.freqs;
    tmpWAVE = interp1(F,tmpWAVE,Fout);
  end
  
  WAVE(:,:,i) = tmpWAVE;
end

out = MatTSA.tfDecomp(WAVE,'decompType','wavelet',...
                           'dataType','complex',...
                           'tVals',tseries.tVals,...
                           'fVals',Fout,...
                           'chanLabels',tseries.chanLabels(dataChans));

                         
waveParams.SCALE = SCALE;
waveParams.COI   = COI;
out.decompParams = varargin;
out.decompParams = [out.decompParams {waveParams}];


end

%% Multi-taper time-frequency decomposition using pmtm
function out = multitaper(tseries,varargin)
% Compute a time-frequency decomposition using Thompson's Multitaper
%
%  out = multitaper(tseries,varargin)
%
% Inputs
% ------
%   tseries : A MatTSA.timeseries object
%
% Param-Value Inputs
% ------------------
%    windowSize : Size of the time window (DEFAULT: 1024)
%            nw : Multitaper Parameter (DEFAULT: 3)
%     FFTLength : Length of FFT (DEFAULT: 1024)
%         freqs : Frequencies to compute decomposition at.
%       nOutput : Number of times to output at
%                   If nOutput<=1: Treated as a fraction of the total
%                                   samples
%                   Otherwise: Treated as an explicit number of samples
%
% REQUIRES THE SIGNAL PROCESSING TOOLBOX
%
% 


%% Input Parsing
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('tseries',@(x) isa(x,'MatTSA.timeseries'));
p.addParameter('windowSize',1024,@(x) isscalar(x)&&isnumeric(x));
p.addParameter(        'nw',   3,@(x) isscalar(x)&&isnumeric(x));
p.addParameter( 'FFTlength',1024,@(x) isscalar(x)&&isnumeric(x));
p.addParameter('freqs',[]);
p.addParameter(   'nOutput',   1,@(x) isscalar(x)&&isnumeric(x));
p.parse(tseries,varargin{:});

%% Execute

% Determine number of outputs
winSize = 2^nextpow2(p.Results.windowSize);
if p.Results.nOutput<=1
  % As a fraction of input length
  nOutput = p.Results.nOutput*size(tseries,1);
else
  % As an explicit number of samples
  nOutput = p.Results.nOutput;
end;
nOutput = floor(nOutput);

% Find window centers
winCenters = linspace(winSize/2,size(tseries,1)-winSize/2,nOutput);
winCenters = ceil(winCenters);

% Eliminate Repeats
winCenters = unique(winCenters);

times = tseries.tRange(1) + (winCenters-1)/tseries.sampleRate;

% Get indices for each window
indices = repmat([-winSize/2+1:winSize/2]',[1 length(winCenters)]);
indices = indices + repmat(winCenters,[size(indices,1) 1]);

% Select frequencies, either using the length of the FFT, or a specific
% list of desired frequencies.
f = p.Results.FFTlength;
if ~isempty(p.Results.freqs)
  f = p.Results.freqs;
end;

% Rearrange data matrix

doChans = tseries.getChannelsByType('data');

for idxChan = 1:numel(doChans)
  disp(['Computing decomposition for channel #' num2str(idxChan)]);
  tseriesData = tseries.data(:,idxChan);
  tseriesData = tseriesData(indices);
    
  % Compute multi-taper
  [pxx,fx] = pmtm(tseriesData,p.Results.nw,f,tseries.sampleRate);
  
  pxxOut(:,:,idxChan) = pxx;
  
end;

out = MatTSA.tfDecomp('multitaper',pxxOut,times,fx,tseries.chanLabels(doChans));
out.params.windowSize = winSize;
out.params.nOutput = nOutput;
out.params.FFTlength = p.Results.FFTlength;
out.params.nw = p.Results.nw;

%out.type = 'multitaper';
%out.tfX = pxx;
%out.fx = fx;
%out.tx = times;

end

%% Fast Fourier Transform based decomposition
function out = fft(tseries,varargin)
% Compute time-frequency decomposition using the FFT
%
% This functionality is likely unnecessary, as the spectrogram method above
% uses the FFT.
%
error('NOT COMPLETE');

taperTypes = {'hanning' 'hamming' 'blackmanharris' 'none'};
p = inputParser;
p.addRequired('tseries',@(x) isa(x,'MatTSA.timeseries'));
p.addParameter('windowSize',[],isscalar(x)&&isnumeric(x));
p.addParameter('ffttaper','hanning',@(x) ismember(x,taperTypes));
p.addParameter('fftlength',1024,isscalar(x)&&isnumeric(x));
p.parse(tseries,varargin{:});

windowSize = p.Results.windowSize;
if isempty(windowSize)
  windowSize = size(tseries,1);
end;

% Get output frequencies
nFreqs = windowSize/2;

end

%% Use EEGLab's timefreq() function to compute the decomposition
function out = eeglab(tseries,varargin)
%% Use EEGLab's timefreq function.
%
% MAY NOT SUPPORT ALL FUNCTIONALITY IN timefreq()
%
% THIS IS LIKELY DEPRECATED BECAUSE FUNCTIONALITY IS DUPLICATED BY THE
% SPECTROGRAM FUNCTION ABOVE
%

p = inputParser;
p.addRequired('tseries',@(x) isa(x,'MatTSA.timeseries'));
p.addParameter('nOutput', 1000, @(x) isscalar(x)&&isnumeric(x));
p.addParameter('windowSize',1024, @(x) isscalar(x)&&isnumeric(x));
p.addParameter('tlimits',[],@(x) (isnumeric(x)&&isvector(x)&&(numel(x)==2)) );
p.addParameter('timesout',[]);
p.addParameter('detrend','off',@ischar);
p.addParameter('type','phasecoher',@ischar);
p.addParameter('cycles',0, @(x) @(x) (isscalar(x)&&isnumeric(x))||...
                                      (isnumeric(x)&&isvector(x)&&(numel(x)==2)));
p.addParameter('verbose','on',@ischar);
p.addParameter('padratio',1,@(x) isscalar(x)&&isnumeric(x));
p.addParameter('freqs',[0 50],@(x) (isnumeric(x)&&isvector(x)&&(numel(x)==2)) );
p.addParameter('freqscale','linear',@ischar);
p.addParameter('nfreqs',[]);
p.addParameter('timeStretchMarks',[]);
p.addParameter('timeStretchRefs',[]);
p.parse(tseries,varargin{:});


g.winsize = 2^nextpow2(p.Results.windowSize);
if p.Results.nOutput<=1
  % As a fraction of input length
  nOutput = p.Results.nOutput*size(tseries,1);
else
  % As an explicit number of samples
  nOutput = p.Results.nOutput;
end;
nOutput = floor(nOutput);

tmioutopt = { 'ntimesout' nOutput };
g.srate = tseries.sampleRate;

if isempty(p.Results.tlimits)
  g.tlimits = [1 size(tseries,1)];
else
  g.tlimits = p.Results.tlimits;
end;

g.detrend   =  p.Results.detrend;
g.type      = p.Results.type;
g.tlimits   = 1000*tseries.tRange; % EEGLab wants things in milliseconds
g.timesout  = 1000*p.Results.timesout;
g.cycles    = p.Results.cycles;
g.verbose   = p.Results.verbose;
g.padratio  = p.Results.padratio;
g.freqs     = p.Results.freqs;
g.freqscale = p.Results.freqscale;
g.nfreqs    = p.Results.nfreqs;
g.timeStretchMarks = p.Results.timeStretchMarks;
g.timeStretchRefs  = p.Results.timeStretchRefs;
timefreqopts = cell(0);

[alltfX freqs timesout R] = timefreq(tseries.data(:,1)', g.srate, tmioutopt{:}, ...
  'winsize', g.winsize, 'tlimits', g.tlimits, 'timesout', g.timesout', 'detrend', g.detrend, ...
  'itctype', g.type, 'wavelet', g.cycles, 'verbose', g.verbose, ...
  'padratio', g.padratio, 'freqs', g.freqs, 'freqscale', g.freqscale, ...
  'nfreqs', g.nfreqs, 'timestretch', {g.timeStretchMarks', g.timeStretchRefs}, timefreqopts{:});

out = MatTSA.tfDecomp('eeglab',alltfX,timesout/tseries.sampleRate,freqs);

%out.type = 'eeglab';
%out.tfX = alltfX;
%out.fx = freqs;
%out.tx = timesout/tseries.sampleRate;

end
