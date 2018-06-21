function pOut = plotWithDecomp(tseriesIn,varargin)
% Plot a MatTSA.timeseries object and an associated decomposition
%
%

%% Input Parsing
p = inputParser;
p.KeepUnmatched = true;
p.addParameter('Parent',[]);
p.addParameter('decompType','wavelet');
p.addParameter('marks',[]);
p.addParameter('imgRange',[]);
p.addParameter('showChan',1);
p.addParameter('showBand',[]);
p.addParameter('showTimes',[]);
p.addParameter('position',[10 10 2550 950]);
p.addParameter('cmap',guiTools.widget.alphacolor);
p.addParameter('logImg',false);
parse(p,varargin{:});

if isempty(p.Results.Parent)
  % Default to opening a new figure
  parent = figure;
else
  parent = p.Results.Parent;
end;

marks = p.Results.marks;

% Position Figure
f1 = figure(parent); clf;
f1.Position = p.Results.position;

%% Display tfDecomps

% Set channels to display
showChan = p.Results.showChan;
if ~iscell(showChan), showChan = {showChan}; end;
nChan = numel(showChan);
if nChan==0, showChan = {[]}; nChan = 1; end;

% Display The Requested tfDecomps
unitHeight = 1/(nChan+1);
cmap = p.Results.cmap;
for i = 1:nChan
  p1(i) = tseriesIn.decomposition.(p.Results.decompType).show(...
    'Parent',f1,...
    'showBand',p.Results.showBand, ...
    'showTimes',p.Results.showTimes,...
    'showChan',showChan{i},...
    'colormap',cmap,...
    'logImg',p.Results.logImg,...
    'range',p.Results.imgRange,...
    'units','normalized',...
    'position',[0.001 (nChan-i+1)*unitHeight 0.999 0.999*unitHeight]);
  drawnow;
end;

drawnow;

%% Plot the timeseries itself
p2 = tseriesIn.plot('Parent',f1,'units','normalized','position',[0.001 0.001*unitHeight 0.999 0.999*unitHeight]);

if ~isempty(marks)
  a = p2.toggleplot.axes;
  axes(a);
  xVal = marks.startOffset/1000;
  plot([xVal(:) xVal(:)],a.YLim,'r','linewidth',2);
  xVal =(marks.startOffset+marks.durations*1000)/1000;
  plot([xVal(:) xVal(:)],a.YLim,'r','linewidth',2);
end;

pOut.p_tf = p1;
pOut.p_eeg = p2;
pOut.listener = addlistener(p2.toggleplot,'updatedOut',@(h,evt) updateDecomp);

  function updateDecomp
    % Function to update displayed time in each tfDecomp
    for k = 1:numel(p1)
      p1(k).showTimes = p2.toggleplot.xrange;
    end;
  end

end