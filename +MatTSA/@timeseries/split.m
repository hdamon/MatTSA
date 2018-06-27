function plotOut = split(tseries,varargin)
% Split Data Plot With Labels
%
% THIS HELP IS OUT OF DATE AND NEEDS TO BE REWRITTEN
%
% Inputs:
%   data (required): (nSamples X nChannels data matrix)
%   ax   (optional): (Axis to plot to)
%
% Param-Value Inputs:
%   'tVals'    : Values along X-axis (otherwise assumed to be 1:N)
%   'dataRange'   : Range of values to display on Y-Axis
%   'chandisp' :
%   'timedisp' :
%   'sampdisp' :
%   'scale'    : Scaling factor
%   'chanLabels'   : Labels for each channel
%   'plotAll'  : Flag to force plotting of all timepoints. When set to 
%                 false, sequences with more than 10k points will be 
%                 subsampled to improve plotting speed.
%
% Written By: Damon Hyde
% Last Edited; May 24, 2016
% Part of the cnlEEG Project
%


%% Input Parsing
p = inputParser;
p.addRequired('tseries',@(x) isa(x,'MatTSA.timeseries'));
p.addOptional('ax',[],@(x) ishghandle(x)&&strcmpi(get(x,'type'),'axes'));
p.addParamValue('tRange',tseries.tRange,@(x) validateattributes(x,{'numeric'},{'vector','numel',2}));
p.addParamValue('dataRange',tseries.dataRange,@(x) validateattributes(x,{'numeric'},{'vector','numel',2}));
p.addParamValue('chandisp',[],@(x) validateattributes(x,{'numeric'},{'vector'}));
p.addParamValue('timedisp',[],@(x) validateattributes(x,{'numeric'},{'vector','numel',2}));
p.addParamValue('sampdisp',[],@(x) validateattributes(x,{'numeric'},{'vector','numel',2}));
p.addParamValue('scale',1,@(x) isnumeric(x)&&numel(x)==1);
p.addParamValue('plotAll',false,@(x) islogical(x));
p.parse(tseries,varargin{:});

ax = p.Results.ax;
tRange = p.Results.tRange;
dataRange = p.Results.dataRange;
scale = p.Results.scale;
tVals = tseries.tVals;
chanLabels = tseries.chanLabels;

%% If no axis provided, open a new figure with Axes
if isempty(ax), figure; ax = axes; end;
axes(ax);

if ~isempty(p.Results.timedisp)
  error('Not yet implemented');
end;

% Get the range of samples to display
if ~isempty(p.Results.sampdisp)
  sampRange = p.Results.sampdisp;
else
  % Default to displaying everything
  sampRange = [1 size(tseries,1)];
end;

if ~isempty(p.Results.chandisp)
  chanDisp = p.Results.chandisp;
else
  chanDisp = ':';
end;

%% For long time series, only render a subset of timepoints      
useIdx = round(linspace(sampRange(1),sampRange(2),10000));
useIdx = unique(useIdx);

tVals = tVals(useIdx);
tRange = [tVals(1) tVals(end)];
chanLabels = chanLabels(chanDisp);

% Get data range and scale
%delta = dataRange(2)-dataRange(1);
delta = max(abs(dataRange));

data = tseries.getPlotData;
data = data(useIdx,chanDisp);

dataChan = tseries.isChannelType('data');
dataChan = dataChan(chanDisp);
data(:,dataChan) = scale * ( data(:,dataChan)./delta );

auxChan = tseries.isChannelType('aux');
auxChan = auxChan(chanDisp);
data(:,auxChan) = scale * ( data(:,auxChan)./delta );

data(:,~dataChan&~auxChan) = data(:,~dataChan&~auxChan)./delta;

% Plot things.
hold on;
for i = 1:size(data,2)
  offset = size(data,2) - (i -1);
  ax.NextPlot = 'add';
  color = 'k'; % Default to black
  if ~dataChan(i)
    color = 'b';  
  end;
  plotOut(i) = plot(tVals,data(:,i)+offset,color,...
                      'ButtonDownFcn',get(ax,'ButtonDownFcn'));
 % set(ax,'ButtonDownFcn',get(plotOut(i),'ButtonDownFcn'));
end;
ax.XLim = tRange;
ax.YLim = [0 size(data,2)+1];

ticks = 1:size(data,2);
if isempty(chanLabels), chanLabels = 1:size(data,2); end
set(ax,'YTick',ticks);
if exist('flip')
  set(ax,'YTickLabel',flip(chanLabels));
else
  % Compatibility w/ earlier matlab versions.
  set(ax,'YTickLabel',flipdim(chanLabels,2));
end
end

