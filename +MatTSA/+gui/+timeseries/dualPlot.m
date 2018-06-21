classdef dualPlot < crlEEG.gui.uipanel
  % Tool for exploring time series timeseries
  %
  % classdef dualPlot < uitools.baseobjs.gui
  %
  % Written By: Damon Hyde
  % Last Edited: June 10, 2016
  % Part of the cnlEEG Project
  %
  
  
  properties (Dependent=true)
    currIdx
    externalSync
  end
    
  properties (Dependent=true, Hidden=true)
    selectedX    % Deprecated
  end
  
  properties (SetAccess=protected)
    miniplot
    toggleplot
    playcontrols
    chanselect
    chanselectbutton
    displayChannels
    vLine
  end
  
  methods
    
    function obj = dualPlot(timeseries,varargin)      
      p = inputParser;
      p.KeepUnmatched = true;
      p.addRequired('timeseries',@(x) isa(x,'crlEEG.type.timeseries'));
      p.addParameter('Parent',[]);
      p.addParamValue('title','TITLE', @(x) ischar(x));
      parse(p,timeseries,varargin{:});           
                
      %% Initialize Base Object
      
      parent = p.Results.Parent;
      if isempty(parent), parent = figure; end;
      
      obj = obj@crlEEG.gui.uipanel('parent',parent,p.Unmatched); %...     
      crlEEG.gui.util.setMinFigSize(gcf,obj.Position(1:2),obj.Position(3:4),5);
                 
      % Set up the channel selection object
      obj.chanselect       = crlEEG.gui.util.selectChannelsFromTimeseries;
      obj.chanselectbutton = obj.chanselect.editButton;
             
      %% Initialize Mini Plot
      obj.miniplot = crlEEG.gui.timeseries.interface.windowPlot(...
        obj.chanselect.output,...
        'parent',obj.panel,...       
        'BorderType','none');
                               
      %% Initialize Toggle Plot
      obj.toggleplot = crlEEG.gui.timeseries.interface.togglePlot(...
                                obj.chanselect.output,...
                                'parent',obj.panel,...
                                'BorderType','none');
            
      %% Initialize Time Controls
      obj.playcontrols = crlEEG.gui.widget.timeplay(...
        'parent',obj.panel,...
        'units','pixels',...
        'position',[325 5 400 30],...                
        'range',[1 size(p.Results.timeseries,1)]);
      obj.playcontrols.BorderType = 'none';
                    
      %% Add listeners 
      obj.setListeners;
                 
      uistack(obj.chanselectbutton,'top');
            
      %obj.size  = p.Results.size;
      %obj.Units = p.Results.units;
      %drawnow;     
      %pause(0.1);      
      %setUnmatched(obj,p.Unmatched);
      
      %crlEEG.gui.util.setMinFigSize(gcf,obj);      
      %set(obj,'Units','normalized');
      
      obj.ResizeFcn = @(h,evt) obj.resizeInternals;
      obj.resizeInternals;
      
      % Copy without decompositions to speed things up.
      tmp = p.Results.timeseries.copy;
      tmp.decomposition = [];
      
      obj.chanselect.input = tmp;
      obj.miniplot.windowEnd = size(p.Results.timeseries,1);
      
    end
    
    function setListeners(obj)
      % Update the toggleplot whenever the miniplot updates
      obj.listenTo{end+1} = ...
        addlistener(obj.miniplot,'updatedOut',@(h,evt) obj.updateToggle);
      
      % Update the miniplot whenever the selected channels are changed.
      obj.listenTo{end+1} = ...
        addlistener(obj.chanselect,'updatedOut',@(h,evt) obj.updateMini);
      
      % Update the current index whenever the automated player updates
      obj.listenTo{end+1} = ...
        addlistener(obj.playcontrols,'updatedOut',@(h,evt)obj.updatedIdx);      
      
      % Update the displayed line whenever the toggleplot updates
      obj.listenTo{end+1} = ...
        addlistener(obj.toggleplot,'updatedOut',@(h,evt) obj.updateLine);    
      
      set(obj.toggleplot.axes,'ButtonDownFcn', @obj.captureMouseClick);
      
      set(get(obj.panel,'Parent'),'KeyPressFcn',@obj.keyPress);      
      
    end
    
    function resizeInternals(obj)      
      % Position the internal components of the dualPlot
      %

      drawnow; % Clear drawing buffer.
      
      % Get current location in pixels
      currUnits = obj.Units;
      obj.Units = 'pixels';
      pixPos = obj.Position;
            
      % Miniplot Position
      miniPlotHeight = max([75/pixPos(4) min([100/pixPos(4) 0.3])]);
      obj.miniplot.Units = 'normalized';
      obj.miniplot.Position = [0 0 1 miniPlotHeight];
      
      % Toggleplot
      togglePlotHeight = max([0.001 1-miniPlotHeight]);
      obj.toggleplot.Units = 'normalized';
      obj.toggleplot.Position = [0 miniPlotHeight 1 togglePlotHeight];
      
      % PlayControls
      
      % This updates the location at which the channelselect GUI will open.
      figPos = getpixelposition(ancestor(obj,'figure'));
      newPos = [figPos(1)-110 figPos(2) 110 figPos(4)];
      obj.chanselect.setPos = newPos;
      
      % Set the correct position for the "Select Channels Button"
      tpPos  = getpixelposition(obj.toggleplot.panel);
      btnPos = getpixelposition(obj.toggleplot.toggleBtn);      
      obj.chanselectbutton.Position = [btnPos(1)+tpPos(1)+105 btnPos(2)+tpPos(2)-1 120 20];
      
      obj.Units = currUnits;      
    end
    
    function set.externalSync(obj,val)
      if val
        obj.playcontrols.externalSync = true;
        obj.listenTo{2}.Recursive = true;
      else
        obj.playcontrols.externalSync = false;
        obj.listenTo{2}.Recursive = false;
      end
    end
    
    function out = get.externalSync(obj)
      out = obj.playcontrols.externalSync;
    end;
    
    function nextStep(obj)
      obj.playcontrols.nextStep;
    end;
    
    function out = get.selectedX(obj)
      warning('dualPlot.selectedX is a deprecated property. Use dualPlot.currIdx instead');
      out = obj.currIdx;
    end
    
    function out = get.currIdx(obj)
      out = [];
      if ~isempty(obj.playcontrols)
      out = obj.playcontrols.currIdx;      
      end;
    end;
    
    function updatedIdx(obj)
      % Callback when obj.playcontrols.currIdx updates
      %
      % Updates the position of the displayed vertical line, and shifts the
      % display window to center that line when it moves out of the
      % currently displayed range.
      %
      shiftLeft  = obj.playcontrols.currIdx<obj.miniplot.windowStart;
      shiftRight = obj.playcontrols.currIdx>obj.miniplot.windowEnd;
      if ~(shiftLeft||shiftRight)
        obj.updateLine;
      elseif shiftLeft
        obj.miniplot.shiftWindow(-5);
      elseif shiftRight
        obj.miniplot.shiftWindow(+5);
      end            
      notify(obj,'updatedOut');
      
    end
    
    function captureMouseClick(obj,h,varargin)
      % ButtonDwnFcn callback for the main data plot axes
      %
      %
      
      a = ancestor(obj.toggleplot.axes,'Figure');
      selType = a.SelectionType;
      switch lower(selType)
        case 'normal'
          obj.pickXVal;
        case 'alt'         
          c = uicontextmenu('Visible','on');
          obj.toggleplot.axes.UIContextMenu = c;
          m1 = uimenu(c,'Label','Mark Event');
          m2 = uimenu(c,'Label','bar');
          m3 = uimenu(c,'Label','baz');
          m(4) = uimenu(c,'Label','Cancel');
        case 'open'
          disp('doubleclick');
        otherwise
          keyboard;
          error('Unknown selection type');
      end      
    end
    
    function pickXVal(obj)
      % ButtonDwnFcn callback for the main data plot axes
      %
      
      pos = get(obj.toggleplot.axes,'CurrentPoint');
      %disp(['Current Position: ' num2str(pos(1))]);
      xvals = obj.miniplot.windowSeries.xvals;
      idx = find(abs(xvals-pos(1))==min(abs(xvals-pos(1))));
      idx = idx(1);
      range = obj.miniplot.windowStart:obj.miniplot.windowEnd;
      %disp(['Current Idx: ' num2str(range(idx))]);
      obj.playcontrols.currIdx = range(idx);         
    end
    
    
    function updateLine(obj)
      % Update (plot) the vertical line denoting the currently selected
      % timepoint.
      range = obj.miniplot.windowStart:obj.miniplot.windowEnd;
      if ~isempty(obj.vLine)&&ishandle(obj.vLine)        
        delete(obj.vLine);        
      end;
      if ~isempty(obj.playcontrols.currIdx)
        if ismember(obj.playcontrols.currIdx,range)
          xVal = obj.miniplot.timeseries.xvals(obj.playcontrols.currIdx);
          yRange = get(obj.toggleplot.axes,'YLim');
          axes(obj.toggleplot.axes); hold on;
          tmp = get(obj.toggleplot.axes,'ButtonDownFcn');
          obj.vLine = plot([xVal xVal],yRange,'r','linewidth',2,...
            'linestyle','--','ButtonDownFcn',tmp);
          set(obj.toggleplot.axes,'ButtonDownFcn',tmp);
          hold off;
            
        end
      end      
    end
    
    function updateMini(obj)      
      
      % Callback to update miniplot input when channel selection is changed
      obj.miniplot.timeseries = obj.chanselect.output;      
    end
    
    function updateToggle(obj)      
      % Callback to update toggle plot when miniplot is changed
      obj.toggleplot.timeseries = obj.miniplot.windowSeries;
      obj.toggleplot.yrange = obj.chanselect.output.yrange;
    end;
    
    function updateImage(obj)
      % Update the data set associated with the toggleplot.
      %
      error('DEPRECATED');
      % Get windowed timeseries and update toggleplot
      timeseries = obj.miniplot.windowSeries;                 
      obj.toggleplot.timeseries = timeseries;     
      
      % Update toggleplot with range of full dataset
      obj.toggleplot.yrange = obj.chanselect.output.yrange;                      
    end
    
    function keyPress(obj,h,evt)
      % function keyPress(obj,h,evt)
      %
      % Callback for handling key controls in the main plot figure;
      %
            
      % Get modifier string
      if isempty(evt.Modifier)
       str = '';
      else
       str = evt.Modifier{1};
      end;
      
      switch evt.Key
        case 'downarrow'
          obj.toggleplot.scale = obj.toggleplot.scale*0.9;          
        case 'uparrow'
          obj.toggleplot.scale = obj.toggleplot.scale*1.1;          
        case 'leftarrow'
          switch str
            case 'shift'
              tic
              obj.miniplot.shiftWindow(-10);
            otherwise          
              tic
              obj.miniplot.shiftWindow(-1);
          end;
        case 'rightarrow'
          switch str
            case 'shift'
              tic
              obj.miniplot.shiftWindow(+10);
            otherwise
              tic
              obj.miniplot.shiftWindow(+1);
          end;
        case 'pageup'
          
          switch str
            case 'shift'           
              % Ensure the number of display channels is always increased
              % by at least 1
              newdisplay = round(1.1*obj.toggleplot.nDisplay);
              if newdisplay ==obj.toggleplot.nDisplay,
                newdisplay = newdisplay + 1;
              end;
              obj.toggleplot.nDisplay = newdisplay;
            otherwise              
              obj.toggleplot.shiftDisplayed(-1);
          end;
        case 'pagedown'
          if isempty(evt.Modifier)
            str = '';
          else
            str = evt.Modifier{1};
          end;
            
          switch str
            case 'shift'     
              % Ensure the number of displayed channels is always decreased
              % by at least 1.
              newdisplay = round(0.9*obj.toggleplot.nDisplay);
              if newdisplay==obj.toggleplot.nDisplay
                newdisplay = newdisplay - 1;
              end
              obj.toggleplot.nDisplay = newdisplay;
            otherwise
             obj.toggleplot.shiftDisplayed(1);
          end;
      end
    end    
    
    
  end
  
  %%  PROTECTED METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Access=protected)
  end
  
  %%  STATIC PROTECTED METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Static=true,Access=protected)
    function p = parseInputs(varargin)

    end
  end
  
end
    