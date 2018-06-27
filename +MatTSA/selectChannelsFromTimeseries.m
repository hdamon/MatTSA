classdef selectChannelsFromTimeseries < handle
%% THIS OBJECT SHOULD BE MOVED OUTSIDE THIS CLASS
%%
  % Select a subset of channels from a MatTSA.timeseries object
  %
  % obj = crlBase.gui.util.selectChannelsFromTimeseries(timeseries)
  %
  % Input
  % -----
  %   timeseries : MatTSA.timeseries object
  %
  % Properties
  % ----------
  %    currChannels : List of currently selected channels
  %    output : MatTSA.timeseries object containing only the 
  %               selected channels.
  %
  % Part of the crlBase Project
  % 2009-2018
  %
  
  properties    
    currChannels
  end
  
  properties (Dependent = true)
    input
    output
  end
  
  properties (Hidden =true, Dependent = true)
    setPos
  end
  
  properties (Access=private)
    gui
    input_
    setPos_ = [2000 100 110 550]; % In Pixels
    outputInternal
    origChannels
  end
  
  events
    updatedOut
  end
  
  methods
    
    function obj = selectChannelsFromTimeseries(timeseries)      
      if nargin>0       
       obj.input = timeseries;
       obj.currChannels = obj.input.chanLabels;
      end
    end
    
    function btnOut = editButton(obj,varargin)
      
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('Position',[1 1 120 20]);
      p.parse(varargin{:});
      
      pos = p.Results.Position;
      pos(3:4) = [50 20];
      
      btnOut = uicontrol('Style','pushbutton',...
                         'String','Select Channels',...
                         'Units','pixels',...
                         'BusyAction','cancel',...
                         'Position',pos,...
                         'CallBack',@(h,evt) obj.editChannels,...
                         p.Unmatched);
                       
    end
    
    function editChannels(obj)
      % Raise or open a GUI to edit selected channels
      %
      if ~guiTools.util.parentfigure.raise(obj.gui)
        
        f = figure('Units','pixels','Position',obj.setPos);
        obj.gui = uicontrol(f,'Style','listbox',...
          'Units','normalized',...
          'Position',[0.02 0.02 0.96 0.96],...
          'Callback',@(h,evt) obj.changeSel);
        set(obj.gui,'Units','normalized');
        obj.syncGUI;
      end
    end
        
    function set.currChannels(obj,val)
      % Set method for
      % crlBase.gui.util.selectChannelsFromTimeseries.selectedStrings
      %
      % Provides input checking, and notifies output when set value is
      % changed (but not otherwise)
      %
      assert(ischar(val)||iscellstr(val),'Input must be a cell string');
      if ischar(val), val = {val}; end;
      if ~isequal(obj.currChannels,val)
        obj.currChannels = val;
        obj.syncGUI;        
        notify(obj,'updatedOut');
      end
    end    
    
  %% GET/SET METHODS FOR DEPENDENT PROPERTIES
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Get/Set Methods for obj.input
    function out = get.input(obj)
      out = obj.input_;
    end
    
    function set.input(obj,timeseries)
      assert(isa(timeseries,'MatTSA.timeseries'),...
              'Input must be a MatTSA.timeseries object');
      if ~isequal(obj.input_,timeseries)
        obj.input_ = timeseries;
        
        if ~isequal(obj.input_.chanLabels,obj.origChannels)
          % Only update the selected channels if the overall list has
          % changed.
          obj.origChannels = obj.input_.chanLabels;
          obj.currChannels = obj.origChannels;
        else          
          notify(obj,'updatedOut');
        end;        
      end      
    end
    
    %% Set/Get methods for obj.output
    function set.output(~,~)
      error('Output of crlBase.gui.util.selectChannelsFromTimeseries is a dependent property');
    end;
    
    function out = get.output(obj)
      out = obj.input(:,obj.currChannels);   
    end    
    
    %% Get/Set methods for obj.setPos
    function out = get.setPos(obj)
      out = obj.setPos_;
    end;
    
    function set.setPos(obj,val)
      % If there's a GUI open, adjust its position.
      obj.setPos_ = val;
      if ~isempty(obj.gui)&&ishghandle(obj.gui)
        set(ancestor(obj.gui,'figure'),'Position',obj.setPos);
      end
    end    
               
  end
  
  %% PRIVATE METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Access=private)
    
    function syncGUI(obj)
      % Make sure that the selected channels in the GUI match those in the
      % obj.currChannels
      %disp('SyncingGUI')
      if ~isempty(obj.input)
        if ~isempty(obj.gui)&&ishghandle(obj.gui)
          allChan = obj.input.chanLabels;
          if ~isequal(allChan,obj.gui.String)
            obj.gui.String = allChan;
            obj.gui.Value = 1:numel(allChan);
            obj.gui.Max = numel(allChan);
          end
          
          chanInInput = find(ismember(obj.input.chanLabels,obj.currChannels));
          if ~isequal(chanInInput,obj.gui.Value)
            obj.gui.Value = chanInInput;
          end;
        end
      end;
    end
    
    function changeSel(obj)            
      % Callback function when obj.gui selection changes
      obj.currChannels = obj.input.chanLabels(obj.gui.Value);
    end;
    
    function closeGUI(obj)      
      % Callback to close gui figure.
        crlBase.gui.util.parentfigure.close(obj.gui);      
    end
  end
  
end
