classdef topoPlot < crlEEG.gui.uipanel
  
  properties
    timeseries
    timepoints
    X
    Y
  end
  
  properties
    axes
    topo
  end
      
  methods
   
    function obj = topoPlot(timeseries,varargin)
      
       %% Input Parsing
      p = inputParser;
      p.KeepUnmatched = true;
      p.addRequired('timeseries',@(x) isa(x,'crlEEG.type.timeseries'));      
      p.addRequired('timepoints',@(x) isnumeric(x)&&isvector(x));      
      p.addOptional('ax',[],@(x) ishghandle(x)&&strcmpi(get(x,'type'),'axes'));      
      p.addParamValue('headNet',[]);
      p.addParamValue('X',[]);
      p.addParamValue('Y',[]);
      p.addParamValue('yrange',[],@(x) isvector(x)&&(numel(x)==2));
      p.addParamValue('scale',0.5,@(x) isnumeric(x)&&numel(x)==1);
           
      parse(p,timeseries,varargin{:});
      
      obj = obj@crlEEG.gui.uipanel(...
        'units','pixels',...
        'position',[10 10 400 400]);            
      
      if ~isempty(p.Results.headNet)
        [x,y] = p.Results.headNet.projPos;
        idx = p.Results.headNet.electrodes.getNumericIndex(timeseries.label);        
        obj.X = x(idx);
        obj.Y = y(idx);
      else
        obj.X = p.Results.X;
        obj.Y = p.Results.Y;
      end;
      
      obj.timeseries = p.Results.timeseries;
      obj.timepoints = p.Results.timepoints;
                  
      obj.doPlot;                      
      
    end
    
    
    function doPlot(obj)
      
      % Nothing to do if no points selected
      if isempty(obj.timepoints), return; end;
      
      nPlot = numel(obj.timepoints);
      
      % Configure Axes
      if isempty(obj.axes)||numel(obj.axes)~=nPlot
        if ~isempty(obj.axes),delete(obj.axes); end;
        
        for i = 1:nPlot
          obj.axes(i) = axes('Parent',obj.panel);
        end;
      end
      
      % Plot TopoPlots
      for i = 1:nPlot
        obj.topo(i) = crlEEG.type.timeseries.render.topo(...
                                         obj.timeseries,...
                                         obj.timepoints(i),...
                                         'x',obj.X,...
                                         'y',obj.Y,...
                                         'ax',obj.axes(i));
      end;
      
    end
    
  end
  
end