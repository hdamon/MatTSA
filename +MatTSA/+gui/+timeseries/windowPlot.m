classdef windowPlot < guiTools.uipanel
  % Plot for Selecting a Sub-Window of Time Data
  %
  % classdef windowPlot < guiTools.uipanel
  %
  % p = uitools.plots.windowPlot(data,varargin)
  %
  % Inputs:
  %   data : nSamples x nChannels Data Matrix
  %
  % Param-Value Pairs:
  %   tVals : x-Axis Values for Each of the N Samples
  %
  % Events:
  %   updatedOut : Notified when the displayed image is updated through
  %                 a call to updateImage.
  %
  % Written By: Damon Hyde
  % Last Edited: June 8, 2016
  % Part of the cnlEEG Project
  %
  
  properties (Dependent=true)
    timeseries
  end
  
  properties (Hidden = true)
    axes
  end
  
  properties (Dependent = true)
    % Tracking of Window Size/Location
    windowStart % Value stored in windowStart_
    windowEnd   % Value stored in windowEnd_
    windowSize  % Computed from windowStart and windowEnd
    
    windowRange
    windowIdx
    windowSeries
  end
  
  properties (Dependent = true, Hidden = true)
    windowData   % Deprecated. Replacement is obj.windowSeries
    windowXVals  % Deprecated. Replacement is obj.windowSeries.tVals
  end
  
  properties (Hidden=true,SetAccess=protected)
    timeseries_
    windowStart_
    windowEnd_
    
    storedVals
    
    zoomButtons
    zoomScale = 0.1;
    
    shiftButtons
    shiftScale = 0.1;
    
    selButtons
    
    vLine
    
  end
  
  methods
    function obj = windowPlot(timeseries,varargin)
      %% uitools.plots.windowPlot object constructor
      
      %% Input Parsing
      p = inputParser;
      p.KeepUnmatched = true;
      p.addRequired('timeseries',@(x) isempty(x)||isa(x,'MatTSA.timeseries'));
      parse(p,timeseries,varargin{:});
      
      %% Initialize Base Object
      obj = obj@guiTools.uipanel(...
        'units','pixels',...
        'position',[10 10 500 150]);
      obj.ResizeFcn = @(h,evt) obj.resizeInternals;
      
      %% Initialize Axes
      obj.axes = axes('Parent',obj.panel);
      
      %% Initialize Zoom Buttons
      obj.zoomButtons = guiTools.widget.dualbutton(...
        'parent',obj.panel,...
        'leftLabel','-',...
        'rightLabel','+',...
        'title','Zoom');
      obj.listenTo{end+1} = addlistener(obj.zoomButtons,'leftPushed',...
        @(h,evt) obj.adjustZoom(+1));
      obj.listenTo{end+1} = addlistener(obj.zoomButtons,'rightPushed',...
        @(h,evt) obj.adjustZoom(-1));
      obj.zoomButtons.BorderType = 'none';
      
      %% Initialize Shift Buttons
      obj.shiftButtons = guiTools.widget.dualbutton(...
        'parent',obj.panel,...
        'title','Shift');
      obj.listenTo{end+1} = addlistener(obj.shiftButtons,'leftPushed',...
        @(h,evt) obj.shiftWindow(-1));
      obj.listenTo{end+1} = addlistener(obj.shiftButtons,'rightPushed',...
        @(h,evt) obj.shiftWindow(+1));
      obj.shiftButtons.BorderType = 'none';
      
      %% Initialize Select Buttons
      obj.selButtons = guiTools.widget.dualbutton(...
        'parent',obj.panel,...
        'title','Select',...
        'leftLabel','Start',...
        'rightLabel','End');
      obj.listenTo{end+1} = addlistener(obj.selButtons,'leftPushed',...
        @(h,evt) obj.setStart(h,evt));
      obj.listenTo{end+1} = addlistener(obj.selButtons,'rightPushed',...
        @(h,evt) obj.setEnd(h,evt));
      obj.selButtons.BorderType = 'none';
      
      obj.resizeInternals;
      
      %% Adjust figure size
      %uitools.setMinFigSize(gcf,obj.origin,obj.size,5);
      
      %% Initialize stored values
      obj.windowStart_ = 1;
      obj.windowEnd_   = 1;
      
      %% Set Data, tVals, and initial windowPlot size, then draw
      obj.timeseries = p.Results.timeseries;
      obj.windowSize = size(obj.timeseries,1);
      
      setUnmatched(obj,p.Unmatched);
      
    end
                
    function resizeInternals(obj)
      % Callback to resize GUI elements when the main size is changed.
      currUnits = obj.Units;
      obj.Units = 'pixels';
      pixPos = obj.Position;
      
      obj.zoomButtons.Units = 'pixels';
      obj.zoomButtons.Position = [5 5 100 35];
      
      obj.shiftButtons.Units = 'pixels';
      obj.shiftButtons.Position = [110 5 100 35];
      
      obj.selButtons.Units = 'pixels';
      obj.selButtons.Position = [215 5 100 35];
      
      axesPos = [5 41 pixPos(3)-10 pixPos(4)-41];
      if axesPos(3)<1, axesPos(3) = 1; end;
      if axesPos(4)<1, axesPos(4) = 1; end;
      
      set(obj.axes,'Units','Pixels');
      set(obj.axes,'Position',axesPos);
            
      % Set units back to where they were
      obj.Units = currUnits;
    end
    
    
    %% Button Callback Functions
    function shiftWindow(obj,val)
      % Shift windowPlot selection left or right.
      %
      % shiftWindow(obj,val)
      %
      % Inputs
      % ------
      %    obj : MatTSA.timeseries.interface.windowPlot object
      %    val : Amount to shift the window.
      %
      % The distance the window is shifted by is computed as:
      %      shift = val*obj.shiftScale*obj.windowSize
      %
      %   obj.windowSize is the size of the current window
      %   obj.shiftScale defaults to 0.1
      %
      % For val > 0 : Window is shifted to the right.
      %     val < 0 : Window is shifted to the left.
      %
      %
      
      %disp('Shifting Window');
      shift = val*obj.shiftScale*obj.windowSize;
      
      % Always shift by at least one sample
      if floor(shift)==0, shift = sign(shift); end;
      
      boundLeft = (obj.windowStart+shift <= 0);
      boundRght = (obj.windowEnd  +shift > size(obj.timeseries,1));
      
      currSize = obj.windowSize;
      if boundLeft
        newWindowStart = 1;
        newWindowEnd   = obj.windowStart + currSize - 1;
      elseif boundRght
        newWindowStart = size(obj.timeseries,1) - currSize + 1;
        newWindowEnd = size(obj.timeseries,1);
      else
        newWindowStart = obj.windowStart + shift;
        newWindowEnd = obj.windowEnd + shift;
      end;
      
      obj.setWindow(newWindowStart,newWindowEnd);
      
      % Prevents wacky things for large shift values
      %       if shift>0
      %         obj.windowEnd = newWindowEnd;
      %         obj.windowStart = newWindowStart;
      %       else
      %         obj.windowStart = newWindowStart;
      %         obj.windowEnd   = newWindowEnd;
      %       end;
      
    end
    
    function adjustZoom(obj,val)
      obj.windowSize = (1+val*obj.zoomScale)*obj.windowSize;
    end
    
    function nearest = sampleNum(obj,val)
      % Return sample number nearest to requested X value
      %
      % nearest = sampleNum(obj,val)
      %
      nearest = abs(val-obj.timeseries.tVals);
      nearest = find(nearest==min(nearest));
    end;
    
    function setEnd(obj,h,evt)
      % Callback for graphically setting end position
      k = 1;
      while k==1
        k = waitforbuttonpress;
      end
      newPos = get(obj.axes,'CurrentPoint');
      
      obj.windowEnd = obj.sampleNum(newPos(1,1));
    end
    
    function setStart(obj,h,evt)
      % Callback for graphically setting start position
      k = 1;
      while k==1
        k = waitforbuttonpress;
      end
      newPos = get(obj.axes,'CurrentPoint');
      
      obj.windowStart = obj.sampleNum(newPos(1,1));
    end;
    
    
    
    %% Methods for Window Size Setting
    function varargout = getWindow(obj)
      % Return the current window
      if nargout==1
       varargout{1} = [obj.windowStart obj.windowEnd];
      else
        varargout{1} = obj.windowStart;
        varargout{2} = obj.windowEnd;
      end;
    end
    
    function setWindow(obj,windowStart,windowEnd)
      % Set the currently selected range for a windowPlot
      %
      % Usage:
      % ------      
      % setWindow(obj,windowPlotObject)
      % setWindow(obj,windowDefinition)
      % setWindow(obj,windowStart,windowEnd)
      %
      % Inputs
      % ------
      %   obj : MatTSA.timeseries.interface.windowPlot object
      %   windowStart : Index to start selection at
      %   windowEnd   : Index to end selection at
      %
      % windowEnd must be greater than windowStart.
      %
      
      % Get start/end from another windowplot object
      if isa(windowStart,'MatTSA.gui.timeseries.interface.windowPlot')
        [windowStart, windowEnd] = windowStart.getWindow;
      end
      
      % Pass in a vector instead of two values
      if numel(windowStart)==2 && ~exist('windowEnd','var')
        windowEnd = windowStart(2);
        windowStart = windowStart(1);
      end;
      
      % Handle Passing Only a Single Parameter
      if ~exist('windowStart','var')||isempty(windowStart)
        windowStart = obj.windowStart_;
      end
            
      if ~exist('windowEnd','var')||isempty(obj.windowEnd)
        windowEnd = obj.windowEnd_;
      end;
      
      % Make sure they're in range
      windowStart = round(min(size(obj.timeseries,1),max(1,windowStart)));
      windowEnd = round(min(size(obj.timeseries,1),max(1,windowEnd)));
      
      % Update Fields
      updated = false;
      if ( obj.windowEnd_~=windowEnd )
        obj.windowEnd_ = windowEnd;
        obj.updateLine(2,obj.windowEnd_);
        updated = true;
      end
      
      if ( obj.windowStart_~=windowStart )
        obj.windowStart_ = windowStart;
        obj.updateLine(1,obj.windowStart_);
        updated = true;
      end;
      
      if updated, notify(obj,'updatedOut'); end;
      
    end
                    
    %% Methods for Getting the Windowed Data
    function out = get.windowSeries(obj)
      out = obj.timeseries(obj.windowStart:obj.windowEnd,:);
    end;
    
    function out = get.windowData(obj)
      % Get method for deprecated property
      warning(['MatTSA.timeseries.interface.windowPlot.windowData '...
        'is deprecated. Use '...
        'MatTSA.timeseries.interface.windowPlot.windowSeries ' ...
        'instead.']);
      out = obj.windowSeries;
    end
    
    function out = get.windowXVals(obj)
      % Get method for deprecated property
      warning('This functionality is deprecated. use obj.windowSeries.tVals instead');
      out = obj.timeseries.tVals(obj.windowStart:obj.windowEnd);
    end
    
    %% Methods for Updating the Data Plot
    function updateImage(obj)
      % Complete redraw of windowPlot
      %
      
      if isempty(obj.timeseries), return; end;
      %tmpFcn = obj.axes.ButtonDownFcn;
      axes(obj.axes);
      cla;
      
      % For long time series, only render a subset of timepoints
      %       tmpSeries = obj.timeseries;
      %       if ( size(tmpSeries,1) > 10000 )
      %         useIdx = round(linspace(1,size(tmpSeries,1),10000));
      %         useIdx = unique(useIdx);
      %         tmpSeries = tmpSeries(useIdx,:);
      %       end;
      
      % Actually render it.
      obj.timeseries.plot('type','butterfly','ax',obj.axes);
      
      %MatTSA.timeseries.render.butterfly(obj.timeseries,'ax',obj.axes);
      
      % Catch zero-length X axes definitions.
      XLim(1) = obj.timeseries.tVals(1); XLim(2) = obj.timeseries.tVals(end);
      if (XLim(1)==XLim(2)),
        XLim(1) = XLim(1) - 0.1;
        XLim(2) = XLim(2) + 0.1;
      end;
      YLim = get(obj.axes,'YLim');
      
      %set(obj.axes,'XTick',[]);
      set(obj.axes,'YTick',[]);
      
      XTick = get(obj.axes,'XTick');
      Label = get(obj.axes,'XTickLabel');
      
      axis(obj.axes);
      for i = 1:numel(XTick)
        text(XTick(i), YLim(1) + 0.1*(YLim(2)-YLim(1)),Label{i});
        obj.axes.XTickLabel{i} = '';
      end;
      %set(obj.axes,'XTickLabel',[]);
      
      axis([XLim(1) XLim(2) YLim(1) YLim(2)]);
      
      obj.drawLines;
    end
    
    
    function updateLine(obj,linenum,idx)
      % Method to update individual vertical lines in a windowPlot
      XVal = obj.timeseries.tVals(idx);
      if numel(obj.vLine)==2
        % Can only update if we've already plotted both.
        set(obj.vLine{linenum},'XData',XVal*[1 1])
      end;
    end;
    
    
    function drawLines(obj)
      % Draw the vertical lines denoting the boundary of the selected
      % region
      %
      
      axes(obj.axes);
      hold on;
      XVal = obj.timeseries.tVals([obj.windowStart obj.windowEnd]);
      YLim = get(obj.axes,'YLim');
      
      obj.vLine{1} = plot(XVal(1)*[1 1],YLim,'r');
      set(obj.vLine{1},'linewidth',2);
      obj.vLine{2} = plot(XVal(2)*[1 1],YLim,'r');
      set(obj.vLine{2},'linewidth',2);
      hold off;
    end
    
    %%  Get/Set Methods for Dependent Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Get/Set Methods for obj.timeseries
    function val = get.timeseries(obj)
      val = obj.timeseries_;
    end;
    
    function set.timeseries(obj,val)
      % Set Method for MatTSA.timeseries.interface.windowPlot
      %
      %
      
      % Input Checking
      if isempty(val), obj.timeseries_ = []; return; end;
      assert(isa(val,'MatTSA.timeseries'),...
              'Input must be a MatTSA.timeseries object');      
            
      if ~isequal(obj.timeseries_,val)
        obj.timeseries_ = val;                
        obj.updateImage;
        notify(obj,'updatedOut');
      end;
    end;
    
    %% Get/Set Methods for obj.windowStart
    function out = get.windowStart(obj)
      out = obj.windowStart_;
    end
    
    function set.windowStart(obj,val)
      
      % Window Always Starts at 1 if data is empty
      if isempty(obj.timeseries), obj.windowStart = 1; return; end;
      
      % Make sure the windowPlot start point is an integer within range
      val = round(min(size(obj.timeseries,1),max(1,val)));
      
      % Update Stored Value
      if val<=obj.windowEnd_&&...
          obj.windowEnd_~=val
        obj.windowStart_ = val;
        obj.updateLine(1,obj.windowStart_);
        notify(obj,'updatedOut');
      end;
      
    end
    
    %% Get/Set Methods for obj.windowEnd
    function out = get.windowEnd(obj)
      out = obj.windowEnd_;
    end;
    
    function set.windowEnd(obj,val)
      % Window always ends at 1 if data is empty
      if isempty(obj.timeseries), obj.windowEnd_ = 1; return; end;
      
      % Make sure windowPlot end point is an integer within range.
      val = round(min(size(obj.timeseries,1),max(1,val)));
      
      if val>=obj.windowStart_&&...
          obj.windowEnd_ ~= val
        obj.windowEnd_ = val;
        obj.updateLine(2,obj.windowEnd_);
        notify(obj,'updatedOut');
      end
    end
    
    %% Get/Set Methods for obj.windowSize
    function out = get.windowSize(obj)
      out = obj.windowEnd_ - obj.windowStart_ + 1;
    end
    
    function set.windowSize(obj,val)
      % Set selection window size
      %
      
      % Do nothing if it's a negative value
      if (val<=0), return; end;
      
      currSize = obj.windowEnd_ - obj.windowStart_ + 1;
      
      delta = val - currSize;
      if abs(delta)<2, delta = sign(delta); end;
      delta = round(delta);
      
      obj.windowStart = obj.windowStart - delta/2;
      obj.windowEnd   = obj.windowEnd + delta/2;
      
    end
    
  end
  
  methods (Static=true,Access=protected)
    function p = parseInputs(varargin)
      
    end
  end
  
end
