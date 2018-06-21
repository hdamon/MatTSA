function varargout = topo(timeseries,timepoint,varargin)
% Topographic plots of timeseries data
%
% function varargout = TOPO(timeseries,timepoint,varargin)
%
% Inputs
% ------
%   timeseries : Input timeseries object
%                 ( Class : crlBase.type.timeseries )
%   timepoint  : Index into timeseries to plot topographic map of.
%   
% Optional Inputs
% ---------------
%   ax : axis to plot to
%
% Param-Value Inputs
% ------------------
%   'headNet' : Headnet for defining the geometry
%                 ( Class : crlBase.head.model.EEG.headNet )
%   'x','y' :  X and Y coordinates to plot the data at. Both must be provided
%                 to use this option.
%
% Part of the crlBase project
% 2009-2017
%


%% Input Parsing
p = inputParser;
p.KeepUnmatched = true;
p.addRequired('timeseries',@(x) isa(x,'crlBase.type.timeseries'));
p.addRequired('timepoint',@(x) isnumeric(x)&&isscalar(x));
p.addOptional('ax',[],@(x) ishghandle(x)&&strcmpi(get(x,'type'),'axes'));
p.addParamValue('headNet',[],@(x) isa(x,'headNet'));
p.addParamValue('x',[],@(x) isnumeric(x)&&isvector(x));
p.addParamValue('y',[],@(x) isnumeric(x)&&isvector(x));
parse(p,timeseries,timepoint,varargin{:});

headNet = p.Results.headNet;

% If no axis provided, open a new figure with Axes
ax = p.Results.ax;
if isempty(ax), figure; ax = axes; end;
axes(ax);

% Get X-Y Positions for plot
if ~isempty(headNet)
  % Plot using a headmap object    
  elecPlot = headNet.plot('axis',ax,'plotchanLabels',timeseries.chanLabels);
  X = elecPlot.scatter.XData;
  Y = elecPlot.scatter.YData;
else  
  assert(~isempty(p.Results.x)&&~isempty(p.Results.y),...
            'Both X and Y values must be provided');
  X = p.Results.x;
  Y = p.Results.y;
  elecPlot.scatter = scatter(X,Y,[],[0 0 0],'filled');
end

tmpUnits = ax.Units;
ax.Units = 'pixels';
minSize = min(ax.Position(3:4));
ax.Units = tmpUnits;

elecPlot.scatter.SizeData = 0.05*minSize;
topoPlot = plotTopo(timeseries.data(timepoint,:),X,Y,p.Unmatched);       
uistack(elecPlot.scatter,'top');

topoOut.figure = gcf;
topoOut.axis = gca;
topoOut.elecPlot = elecPlot;
topoOut.topoPlot = topoPlot;
topoOut.headCartoon = crlBase.gui.util.drawHeadCartoon(gca,'diam',1.2);

if nargout>0
  varargout{1} = topoOut;
end;

end

function varargout = plotTopo(data,x,y,varargin)
  % Actual function for drawing a topographic plot
  %
  % function topoOut = PLOTTOPO(data,x,y,varargin)
  %
  % Inputs
  % ------
  %    data : data to plot
  %     x,y : X-Y locations of data
  %   
  % Param-Value Inputs
  % ------------------
  %  'diam'   : Size of plot square
  %               ( Default : 4 )
  %  'scale'  : Scaling of input data
  %               ( Default : 1 )
  %  'thresh' : Size of circle to plot value inside
  %               ( Default : 1.6 )
  %  'cmap' : Colormap to use 
  %             ( Class : crlBase.gui.widget.alphacolor )
  %
  % Outputs
  % -------
  %   topoOut : Structure containing plotted image and alphacolor map.
  %
  % Part of the crlBase Projects
  % 2009-2017
  %

  % Input Parsing
  p = inputParser;
  p.addParamValue('diam',4);
  p.addParamValue('scale',1);
  p.addParamValue('thresh',1.6);
  p.addParamValue('cmap',[],@(x) isa(x,'crlBase.gui.widget.alphacolor'));
  p.parse(varargin{:});
 
  if isempty(p.Results.cmap)
    cmap = crlBase.gui.widget.alphacolor;
    cmap.range = [min(data(:)) max(data(:))];
  else
    cmap = p.Results.cmap; 
  end;
  topoOut.cmap = cmap;
  
  diam = p.Results.diam;
  scale = p.Results.scale;
  
  x = scale*x; y = scale*y;
  F = scatteredInterpolant(x(:),y(:),data(:));
  %F.ExtrapolationMethod = 'none';
  
  xGrid = linspace(-diam/2,diam/2,500);
  yGrid = linspace(-diam/2,diam/2,500);
  [X Y] = meshgrid(xGrid,yGrid);
  
  D = F(X,Y);
  
  % NaN's outside the circle
  r = sqrt(X.^2+Y.^2); 
  Q = r>p.Results.thresh;
  D(Q) = nan;
  
  [Dimg, Dalpha] = cmap.img2rgb(D);
  
  holdStatus = ishold(gca);
  hold on;
  topoOut.img = image(xGrid,yGrid,Dimg);
  topoOut.img.AlphaData = ones(size(D));
  topoOut.img.AlphaData(isnan(D)) = 0;
  if ~holdStatus, hold off; end;
  
  a = gca; a.YDir = 'normal'; axis off;
  
  if nargout>0
    varargout{1} = topoOut;
  end;
end
