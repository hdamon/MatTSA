classdef togglePlot < guiTools.uipanel
  % An interface for MatTSA.timeseries data 
  %
  % classdef toggle < uitools.cnlUIObj
  %
  % Provides a UI interface for toggling between a split and a butterfly
  % plot. Allows interactive selection of a subset of channels for display.
  %
  % Properties:
  %    dataRange : Initially set to 
  %
  % Settable Properties
  % ---------
  %    timeseries : MatTSA.timeseries object to display
  %    
  %    
  % By default, togglePlot displays all channels in the input timeseries.
  % It can be configured to display only a subset of channels using the
  % following properties:
  %    displayRange : 1x2 array with First/Last channels to display
  %    nDisplay     : Dependent property. When set, adjusts the values
  %                     in displayRange to display the requested number
  %                     of channels.
  % A toggleplot can only display a contiguous block of channels. This is
  % to permit well defined shifting of the display channels through the two
  % buttons on the right hand side.
  %
  % Written By: Damon Hyde
  % Last Edited: June 1, 2016
  % Part of the cnlEEG Project
  %  
  
  properties (SetObservable,AbortSet)
    timeseries    
    scale
    doSplit = false;
  end
  
  properties
    axes
    plot    
  end
  
  properties (Dependent)
    tRange
    dataRange
  end
  
  properties (Hidden=true)
    toggleBtn        
    shiftBtn    
  end
  
  properties (Hidden=true,Dependent)
    displayChannels
    displayRange
    nDisplay
  end
  
  properties (Access=private)
    tRange_
    dataRange_
    chan_    
    numChan
  end
  
  methods
    
    function obj = togglePlot(timeseries,varargin)
   
      %% Input Parsing
      p = inputParser;
      p.KeepUnmatched = true;
      p.addRequired('timeseries',@(x) isempty(x)||isa(x,'MatTSA.timeseries'));      
      p.addOptional('ax',[],@(x) ishghandle(x)&&strcmpi(get(x,'type'),'axes'));      
      p.addParamValue('dataRange',[],@(x) isvector(x)&&(numel(x)==2));
      p.addParamValue('scale',1,@(x) isnumeric(x)&&numel(x)==1);
           
      parse(p,timeseries,varargin{:});
                  
      %% Initialize cnlUIObj
      obj = obj@guiTools.uipanel(...
          'units','pixels',...
          'position',[10 10 600 600]);
      
      obj.ResizeFcn = @(h,evt) obj.resizeToggleplot;
        
      %% Add the Toggle Button
      obj.toggleBtn = uicontrol('Parent',obj.panel,...
        'Style','pushbutton',...
        'String','Split/Butterfly');%,...
        %'Units','pixels',...
        %'Position',[5 5 100 20]);
      set(obj.toggleBtn,'Callback',@(h,evt) obj.toggleSplit);
      %set(obj.toggleBtn,'Units','normalized');
     
                       
      %% Set up Plot Axes
      obj.axes = axes('Parent',obj.panel); %,'Units','Pixels',...
        %'Position',[70 30 300 300]); 
      
      obj.shiftBtn(1) = uicontrol('Parent',obj.panel,...
        'Style','pushbutton',...
        'String',char(8595),...
        ...%'Units','normalized',...
        ...%'Position',[0.96 axisYstart 0.02 axisYsize/2-0.02],...
        'Callback',@(h,evt) obj.shiftDisplayed(1));
      
      obj.shiftBtn(2) = uicontrol('Parent',obj.panel,...
        'Style','pushbutton',...
        'String',char(8593),...
        ...%'Units','normalized',...
        ...%'Position',[0.96 axisYstart+axisYsize/2 0.02 axisYsize/2-0.02],...
        'Callback',@(h,evt) obj.shiftDisplayed(-1));
              
      % Initial Sizing Callback
      obj.Position = [10 10 600 600];
      
      %% Set Property Values
       obj.dataRange = p.Results.dataRange;
      obj.scale  = p.Results.scale;
      obj.timeseries   = p.Results.timeseries;      
           
            
      % Select Channels to Display. Up to 30 is the default.

      
      % Set Desired UIPanel properties
      obj.setUnmatched(p.Unmatched);
      
      %MatTSA.gui.util.setMinFigSize(gcf,obj);
      
      %% Do Initial Display of Plot
      obj.updateImage;
    end
    
    function set.scale(obj,val)
      % Update the image when the scale changes.
      if ~isequal(obj.scale,val)
        obj.scale = val;
        obj.updateImage;
      end;
    end
    
    %% Get/Set dataRange
    function val = get.dataRange(obj)
      if ~isempty(obj.dataRange_)
        val = obj.dataRange_;
      else
        val = obj.timeseries.dataRange;
      end;
    end
    
    function set.dataRange(obj,val)      
      if ~isequal(obj.dataRange_,val)
        obj.dataRange_ = val;        
        obj.updateImage;
      end;       
    end
    
    function val = get.tRange(obj)
      if ~isempty(obj.tRange_)
        val = obj.tRange_;
      else
        val = obj.timeseries.tRange;
      end;
    end
    
    function set.tRange(obj,val)      
      if ~isequal(obj.tRange_,val)
        obj.tRange_ = val;        
        if ishghandle(obj.plot)
          keyboard;
        end;
      end;       
    end
    
    
    %% Set method for internal timeseries
    function set.timeseries(obj,val)
      assert(isa(val,'MatTSA.timeseries'),...
              'Must be a MatTSA.timeseries object');            
     % if ~isequal(obj.timeseries,val)
        % Only update if there's a change.
        obj.timeseries = val;  
        obj.checkChan;
        obj.updateImage;
      %end;
    end
    
    function checkChan(obj)
      if size(obj.timeseries,2)~=obj.numChan
        obj.chan_ = [];
        obj.numChan = size(obj.timeseries,2);
      end;        
    end
    
    %% 
    function resizeToggleplot(obj)
      % Callback to adjust internal sizings of toggleplot uipanel when
      % parent panel is resized.
      
      % Size Definitions
      btnWidth    = 15;
      axesXOffset = 70;
      axesYOffset = 35;
      toggleBtnSize = [ 5 5 100 20];
      
      % Do Everything in Pixels
      currUnits = obj.Units;
      obj.Units = 'pixels';
      pixPos = obj.Position;
      
    %  disp(['EEG panel is of size: ' num2str(pixPos)]);
      
      xSize = pixPos(3);     
      ySize = pixPos(4);
      
      % Toggle Button in Lower Left Corner
      set(obj.toggleBtn,'Units','pixels');
      set(obj.toggleBtn,'Position',toggleBtnSize);                  
      
      % Axes Offset from Toggle Button
      obj.axes.Units = 'pixels';
      axesPos = get(obj.axes,'Position');
      axesPos(1) = axesXOffset;
      axesPos(2) = toggleBtnSize(2)+axesYOffset;
      axesPos(3) = xSize - axesPos(1) - btnWidth - 5; if axesPos(3)<=0, axesPos(3) = 1; end;      
      axesPos(4) = ySize - axesPos(2); if axesPos(4)<=0, axesPos(4) = 1; end;
     % disp(['Setting EEG Axes Position: ' num2str(axesPos)]);
      set(obj.axes,'Position',axesPos);
      
      % Shift Buttons
      btnHeight = (ySize-axesPos(2)-5)/2;      
      if btnHeight <= 0, btnHeight = 1; end;
      set(obj.shiftBtn(1),'Units','pixels');
      set(obj.shiftBtn(1),'Position',[xSize-btnWidth axesPos(2) btnWidth btnHeight]);
      
      set(obj.shiftBtn(2),'Units','pixels');
      set(obj.shiftBtn(2),'Position',[xSize-btnWidth axesPos(2)+btnHeight+5 btnWidth btnHeight]);
      
      %disp(['Panel Position:' num2str(pixPos)]);
      obj.Units = currUnits;
    end
    
    %% Get/Set for obj.nDisplay
    function out = get.nDisplay(obj)
      % Get method for togglePlot.nDisplay
      tmp = obj.displayRange;
      out = tmp(2) - tmp(1) + 1;
    end;
    
    function set.nDisplay(obj,val)
      % Set method for togglePlot.nDisplay
      val = round(val);
      if val==obj.nDisplay
        return;        
      else
        % Change Number of Channels
        currRange = obj.displayRange;
        
        newRange = [currRange(1) currRange(1)+val-1];        
        
        % If we go off the end, add it to the beginning instead
        if newRange(2)>size(obj.timeseries,2)
          extraLen = newRange(2) -size(obj.timeseries,2);
          newRange(2) = size(obj.timeseries,2);
          newRange(1) = newRange(1) - extraLen;
        end
        
        % Can't go before the first index
        if newRange(1)<1, newRange(1) = 1; end;
        
        % Consistency check
        if newRange(2)<newRange(1), newRange(2) = newRange(1); end;
        
        obj.displayRange = newRange;
        %obj.updateImage;
      end;      
    end
    
    %% Get/Set for obj.displayRange
    function set.displayRange(obj,val)
      % Set the range of channels to display in the split toggleplot.
      %
      assert(isnumeric(val)&&(numel(val)==2),...
        'Input must be a numeric vector with two elements');
      
      if ( val(1) < 1 )
        val(1) = 1;
      end;
      
      if ( val(2) > size(obj.timeseries,2) )
        val(2) = size(obj.timeseries,2);
      end;
      
      if ~isequal(val,obj.chan_)
        obj.chan_ = val;
        obj.updateImage;
      end;
    end
    
    function val = get.displayRange(obj)
      if ~isempty(obj.chan_)
        if obj.chan_(2)>size(obj.timeseries,2)
          % This might be unnecessary?
          obj.chan_(2) = size(obj.timeseries,2);
        end;
        val = obj.chan_;
      else
        % Default if it hasn't been set yet.
        nDisp = 30;
        if size(obj.timeseries.data,2)<nDisp, nDisp = size(obj.timeseries.data,2); end;
        val = [1 nDisp];
      end;
    end
    
    function shiftDisplayed(obj,shiftDir,shiftDist)    
      % Shift the window on which channels are displayed
      %
      % SHIFTDISPLAYED(obj,shiftDir,shiftDist)
      %
      % Inputs
      % ------
      %   shiftDir : Direction (Typically +1/-1) to shift display
      %   shiftDist : (Optional) Distance to shift the display as a
      %                 fraction of the number of currently displayed
      %                 channels.  DEFAULT: 0.25

      
      if ~exist('shiftDist','var'), shiftDist = 0.25; end;
      
      shift = shiftDir*round(shiftDist*obj.nDisplay);
      if shift==0, shift = 1; end;
      
      newRange = obj.displayRange +shift;
      
      if newRange(1)<1, 
        newRange = newRange - newRange(1) + 1; 
      end;
      
      if newRange(2)>size(obj.timeseries,2),
        newRange = newRange - (newRange(2)-size(obj.timeseries,2));
      end;
      
      if newRange(2)<newRange(1)
        keyboard;
      end;
      
      obj.displayRange = newRange;            
    end
    
    function set.doSplit(obj,val)
      if obj.doSplit ~= val
         obj.doSplit = val;
         obj.updateImage;
      end;
    end

    %%
    function toggleSplit(obj)
      % Callback for toggle between a split and butterfly plot
      %
      obj.doSplit = ~obj.doSplit;      
    end
      
    %%
    function updateImage(obj)      
      % Update the image displayed in the toggleplot
      %      
      try
        % Clear axis, and make sure the next plot doesn't modify callbacks
        tmpFcn = obj.axes.ButtonDownFcn;
        cla(obj.axes,'reset');
        obj.axes.ButtonDownFcn = tmpFcn;
                
%         tmpSeries = obj.timeseries;
%         if ( size(tmpSeries,1) > 10000 )
%           useIdx = round(linspace(1,size(tmpSeries,1),10000));
%           useIdx = unique(useIdx);
%           tmpSeries = tmpSeries(useIdx,:);
%         end;
        
        % Plot
        if obj.doSplit % Do a split plot
          dispChan = obj.displayRange(1):obj.displayRange(2);          
          
          obj.plot = split(obj.timeseries(:,dispChan),obj.axes,...
            'tRange',obj.tRange,'dataRange',obj.dataRange,'scale',obj.scale);
        else % Just do a butterfly plot
          obj.plot = butterfly(...
                                                obj.timeseries,obj.axes,...
                                                'tRange',obj.tRange,...
                                                'dataRange',obj.dataRange,...
                                                'scale',obj.scale);
        end;
        drawnow;
        notify(obj,'updatedOut');
      catch
        % This is a bit of a stupid thing, because I'm using it to avoid
        % explicitly checking to see whether I can plot or not.
      end;
    end
    
  end
  
end