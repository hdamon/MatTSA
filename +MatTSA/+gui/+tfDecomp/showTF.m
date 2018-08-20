classdef showTF < guiTools.uipanel
% Provides a GUI interface for crlEEG.type.timeFrequencyDecomposition objects
%
% classdef showTF < matlab.ui.container.Panel
% 
% Inputs
% ------
%   tfDecomp : MatTSA.tfDecomp Object to display
%
% Param-Value Inputs
% ------------------
%      'title' : Title for the plot 
%      'range' : Range to use for the colormap
%   'colormap' : guiTools.widget.alphacolor object
%     'logImg' : Flag to enable plotting of a log10-scaled image
%   'showChan' : Cell 
%   'showBand'
%  'showTimes'
%
% Additional parameters will be passed to a guiTools.uipanel object. See
% that help for more information.
%
%
  properties
    tfDecomp   % The decomposition being displayed
    ax         % Handle to the axis the display is in
    cmap       % Colormap object
  end;
  
  properties (Dependent =true)
    showChan
    showBand
    showTimes
    logImg
    imgRange    
  end
  
  % Publically accessible but hidden
  properties (Hidden)
    chanSelect
    editCMap
  end
    
  %% Actual Parameter Values are stored in Protected Properties
  properties (Access=protected)
    showBand_
    showTimes_
    logImg_
    imgRange_ = [];
    listeners_ 
    cbar_
  end
  
  
  %% METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods
    
    function obj = showTF(tfDecomp,varargin)
      
        %% Input Parsing
        p = inputParser;
        p.KeepUnmatched = true;
        p.addRequired('tfDecomp',@(x) isa(x,'MatTSA.tfDecomp'));
        p.addParameter('title','TITLE',@(x) ischar(x));
        p.addParameter('showBand',[]);
        p.addParameter('showTimes',[]);
        p.addParameter('showChan',[]);
        p.addParameter('logImg',false);
        p.addParameter('range',[]);        
        p.addParameter('colormap',guiTools.widget.alphacolor,@(x) isa(x,'guiTools.widget.alphacolor'));
        p.parse(tfDecomp,varargin{:});
                      
        %% Superclass Constructor
        obj = obj@guiTools.uipanel(p.Unmatched); 
                
        %% Display Axes
        obj.ax = axes('parent',obj.panel,'units','normalized');        
        
        %% Interactive Colormap
        obj.cmap = p.Results.colormap;
        if ~isempty(p.Results.range)
          obj.cmap.range = p.Results.range;
        end;
                    
        %% Input Data Handle
        obj.tfDecomp = p.Results.tfDecomp;
        
        %% Channel Selection Object
        obj.chanSelect = uicontrol('Style','popup',...
                                   'String',obj.tfDecomp.chanLabels,...
                                   'Parent',obj.panel,...                                   
                                   'CallBack',@(h,evt) updateImage(obj));
                                 
        obj.editCMap = uicontrol('Style','pushbutton',...
                                  'String','Edit Colormap',...
                                  'Parent',obj.panel,...
                                  'CallBack',@(h,evt) obj.cmap.edit);
                                
        
        if ~isempty(p.Results.showChan)
          if iscellstr(p.Results.showChan)||ischar(p.Results.showChan)
            idx = find(cellfun(@(x) isequal(x,p.Results.showChan),obj.tfDecomp.chanLabels));
          else
            idx = p.Results.showChan;
          end;
          assert(numel(idx)==1,'Only one channel can be displayed at a time');
          obj.chanSelect.Value = idx;
        end
                
        %% Set Internal Variables
        obj.showBand_  = p.Results.showBand;
        obj.showTimes_ = p.Results.showTimes;
        obj.logImg_    = p.Results.logImg;
        
        %% Add Listeners
        obj.listeners_ = addlistener(obj.cmap,'updatedOut',@(h,evt) obj.updateImage);
        
        %% Set Resize CallBack. Then Call it.
        obj.ResizeFcn = @(h,evt) obj.resizeInternals();
        obj.Units = 'normalized';
        obj.resizeInternals;
        
        %% Actually Plot
        obj.updateImage;                
    end      
    
    %% Get/Set Methods for Dependent Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function out = get.showBand(obj)
      out = obj.showBand_;
    end;
    
    function set.showBand(obj,val)
      obj.setProperty('showBand',val);      
    end

    function out = get.showTimes(obj)
      out = obj.showTimes_;
    end
    
    function set.showTimes(obj,val)
      obj.setProperty('showTimes',val);
    end    
    
    function out = get.logImg(obj)
      out = obj.logImg_;
    end;
    
    function set.logImg(obj,val)
      obj.setProperty('logImg',val);      
    end
    
    function out = get.imgRange(obj)
      out = [];
      if ~isempty(obj.cmap)
        out = obj.cmap.range;
      end      
    end
    
    function set.imgRange(obj,val)
      if ~isempty(obj.cmap)
        obj.cmap.range = val;
      end;      
    end
              
    function out = get.showChan(obj)
      out = obj.chanSelect.String{obj.chanSelect.Value};
    end
    
    function set.showChan(obj,val)
      newIdx = find(cellfun(@(x) isequal(x,val),obj.chanSelect.String));
      if numel(newIdx)==0, error('String not found'); end;
      if numel(newIdx)>1, error('Multiple match'); end;
      
      if ~(newIdx==obj.chanSelect.Value)
        obj.chanSelect.Value = newIdx;
        obj.updateImage;
      end
        
    end

    %% GUI Callbacks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function resizeInternals(obj)
      % Callback for resizing 
      %
      %
      currUnits = obj.Units;
      cleanup = onCleanup(@() set(obj,'Units',currUnits));
      
      obj.Units = 'pixels';
      pixPos = obj.Position;    
      obj.chanSelect.Units = 'pixels';
      obj.chanSelect.Position = [2 2 100 30];
            
      obj.editCMap.Units = 'pixels';
      obj.editCMap.Position = [ 105 2 100 30];
      
      obj.ax.Units = 'pixels';
      xSize = max([5 (pixPos(3)-95)]);      
      ySize = max([5 0.95*(pixPos(4)-50)]);      
      obj.ax.Position = [70 50 xSize ySize];                        
    end
    
    function updateImage(obj)
      % Update the displayed image
      %
      
      axes(obj.ax); cla;
      
      % Display the Image
      set(guiTools.util.parentfigure.get(obj),'colormap',obj.cmap.cmap);
      obj.tfDecomp.imagesc('parent',obj.ax,...
                           'showBand',obj.showBand,...
                           'showChan',obj.showChan,...
                           'showTimes',obj.showTimes,...
                           'logImg',obj.logImg,...
                           'colormap',obj.cmap);
                         
      % Generate a Colorbar
      if isempty(obj.cbar_)||~ishghandle(obj.cbar_)
        obj.cbar_ = colorbar('peer',obj.ax,'East');                   
      end;
      obj.cbar_.FontSize = 20;
      obj.cbar_.FontWeight = 'bold';
      
      tvals = obj.cbar_.Ticks;
      tickLabels = strsplit(num2str(obj.cmap.range(1) + tvals*(obj.cmap.range(2)-obj.cmap.range(1))));
      obj.cbar_.TickLabels = tickLabels;
    end
    
  end
  
  %%  PROTECTED METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Access=protected)
    
    function setProperty(obj,prop,val)
      internalProp = [prop '_'];
      if ~isequal(obj.(internalProp),val)
        obj.(internalProp) = val;
        obj.updateImage;
      end;      
    end
    
  end
  
  %%  STATIC PROTECTED METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Static=true,Access=protected)
    function p = parseInputs(varargin)

    end
  end
    
  
end
