function plotOut = butterfly(tseries,varargin)
% Butterfly tseries plot
%
% plotOut = BUTTERFLY(tseries,varargin)
%
% Required Inputs:
%   tseries :  Nsamples x Nchannels Data Matrix
% 
% Optional Inputs:
%   'ax'     : Handle to a matlab axis to display in
%   'tVals'  : Values to display along X-Axis
%   'dataRange' : Range to use for Y-Axis
%
% This object exists mostly to maintain a consistent structure across the
% entire uitools package. It additionally ensures that the resulting plot
% and axes keep the same button down function that the axis had prior to
% doing the plot.
%
% Written By: Damon Hyde
% Last Edited: May 24, 2016
% Part of the cnlEEG Project
%

import crlBase.util.validation.isNumericVector

%% Input Parsing
p = inputParser;
p.addRequired('tseries',@(x) isa(x,'crlBase.type.timeseries'));
p.addOptional('ax',[],@(x) ishghandle(x)&&strcmpi(get(x,'type'),'axes'));
p.addParamValue('tRange',tseries.tRange,@(x) isvector(x)&&(numel(x)==2));
p.addParamValue('dataRange',tseries.dataRange,@(x) isvector(x)&&(numel(x)==2));
p.addParamValue('chandisp',[],@(x) isNumericVector(x));
p.addParamValue('timedisp',[],@(x) isNumericVector(x,2));
p.addParamValue('sampdisp',[],@(x) isNumericVector(x,2));
p.addParamValue('scale',1,@(x) isnumeric(x)&&numel(x)==1);
p.addParamValue('plotAll',false,@(x) islogical(x));
p.parse(tseries,varargin{:});

ax = p.Results.ax;

tVals = tseries.tVals;
tRange = p.Results.tRange;
dataRange = p.Results.dataRange;

%% Shift to appropriate axes, or open a new one.
if isempty(ax), figure; ax = axes; end;
axes(ax);

% If plotting a single timepoint, use X's.
if numel(tVals)==1, plotOpts = 'x';
else                plotOpts = '';  end;

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

%% Plot!
%plotOut = plot(tVals,tseries.data,plotOpts);
tmpData = tseries.getPlotData;
tmpData = tmpData(useIdx,chanDisp);

dataChans = tseries.isChannelType('data');
dataChans = dataChans(chanDisp);

ax.NextPlot = 'add';
if any(dataChans)
  plotOut = plot(tVals(useIdx),tmpData(:,dataChans),['k' plotOpts],'ButtonDownFcn',get(ax,'ButtonDownFcn'));
end

if any(~dataChans)
  ax.NextPlot = 'add';
  tmp = plot(tVals(useIdx),tmpData(:,~dataChans)./p.Results.scale,[plotOpts], ...
            'linewidth',2,'ButtonDownFcn',get(ax,'ButtonDownFcn'));
  if exist('plotOut','var')
    plotOut = [plotOut ; tmp];
  else
    plotOut = tmp;
  end;
end;

% Modify X Limits if plotting a single timepoint.
XLim = p.Results.tRange;
if XLim(1)==XLim(2), XLim = XLim + [-0.1 0.1]; end;
ax.XLim = XLim;

% Make the y limits symmetric.
YLim = max(abs(dataRange));
if YLim(1)==0, YLim = 0.1; end;
YLim = YLim./p.Results.scale;
ax.YLim = [-YLim YLim];


% Set Yticks.
ticks = linspace(0,YLim,6);
ticks = [-flipdim(ticks(2:end),2) ticks];
set(ax,'YTick',ticks);

e = log10(ticks(end));
e = sign(e)*floor(abs(e));
yt = ticks/10^e;

chanLabels = cell(size(yt));
for i = 1:numel(yt)
  chanLabels{i} = sprintf('%1.2f',yt(i));
end

set(ax,'YTickLabel',chanLabels);

text(XLim(1),YLim(1)*0.95,sprintf('\\times 10^{%d}',e));

end

