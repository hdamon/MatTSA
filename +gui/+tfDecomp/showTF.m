classdef showTF < crlEEG.gui.uipanel
%classdef showTF < matlab.ui.container.Panel
% Provides a GUI interface for crlEEG.type.timeFrequencyDecomposition objects
%
% 
  properties
    tfDecomp
    ax
    cmap
  end;
  
  properties (Dependent =true)
    showBand
    showTimes
    logImg
    imgRange
    showChan
  end
   
  properties (Hidden)
    chanSelect
    editCMap
  end
  
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
        p.addRequired('tfDecomp',@(x) isa(x,'timeFrequencyDecomposition'));
        p.addParameter('title','TITLE',@(x) ischar(x));
        p.addParameter('showBand',[]);
        p.addParameter('showTimes',[]);
        p.addParameter('showChan',[]);
        p.addParameter('logImg',false);
        p.addParameter('range',[]);        
        p.addParameter('colormap',crlEEG.gui.widget.alphacolor,@(x) isa(x,'crlEEG.gui.widget.alphacolor'));
        p.parse(tfDecomp,varargin{:});
                      
        % Superclass Constructor
        obj = obj@crlEEG.gui.uipanel(p.Unmatched); 
                
        % Display Axes
        obj.ax = axes('parent',obj.panel,'units','normalized');        
        
        % Interactive Colormap
        obj.cmap = p.Results.colormap;
        if ~isempty(p.Results.range)
          obj.cmap.range = p.Results.range;
        end;
             

        
        % Input Data Handle
        obj.tfDecomp = p.Results.tfDecomp;
        
        % Channel Selection Object
        obj.chanSelect = uicontrol('Style','popup',...
                                   'String',obj.tfDecomp.labels,...
                                   'Parent',obj.panel,...                                   
                                   'CallBack',@(h,evt) updateImage(obj));
                                 
        obj.editCMap = uicontrol('Style','pushbutton',...
                                  'String','Edit Colormap',...
                                  'Parent',obj.panel,...
                                  'CallBack',@(h,evt) obj.cmap.edit);
                                
        
        if ~isempty(p.Results.showChan)
          if iscellstr(p.Results.showChan)||ischar(p.Results.showChan)
            idx = find(cellfun(@(x) isequal(x,p.Results.showChan),obj.tfDecomp.labels));
          else
            idx = p.Results.showChan;
          end;
          assert(numel(idx)==1,'Only one channel can be displayed at a time');
          obj.chanSelect.Value = idx;
        end
                
        % Set Internal Variables
        obj.showBand_  = p.Results.showBand;
        obj.showTimes_ = p.Results.showTimes;
        obj.logImg_    = p.Results.logImg;
        
        % Add Listeners
        obj.listeners_ = addlistener(obj.cmap,'updatedOut',@(h,evt) obj.updateImage);
        
        % Set Resize CallBack. Then Call it.
        obj.ResizeFcn = @(h,evt) obj.resizeInternals();
        obj.Units = 'normalized';
        obj.resizeInternals;
        
        % Actually Plot
        obj.updateImage;                
    end      
    
    function out = get.showBand(obj)
      out = obj.showBand_;
    end;
    
    function set.showBand(obj,val)
      if ~isequal(obj.showBand_,val)
        obj.showBand_ = val;
        obj.updateImage;
      end;
    end

    function out = get.showTimes(obj)
      out = obj.showTimes_;
    end
    
    function set.showTimes(obj,val)
      if ~isequal(obj.showTimes_,val)
        obj.showTimes_ = val;
        obj.updateImage;
      end
    end    
    
    function out = get.logImg(obj)
      out = obj.logImg_;
    end;
    
    function set.logImg(obj,val)
      if ~isequal(obj.logImg_,val)
        obj.logImg_ = val;
        obj.updateImage;
      end;
    end
    
    function out = get.imgRange(obj)
      out = obj.imgRange_;
    end
    
    function set.imgRange(obj,val)
      if ~isequal(obj.imgRange_,val)
        obj.imgRange_ = val;
        obj.updateImage;
      end;
    end
          
    function resizeInternals(obj)
      currUnits = obj.Units;
      cleanup = onCleanup(@() set(obj,'Units',currUnits));
      
      obj.Units = 'pixels';
      pixPos = obj.Position;
    %  disp(['TF Panel is of Size: ' num2str(pixPos)]);
      obj.chanSelect.Units = 'pixels';
      obj.chanSelect.Position = [2 2 100 30];
            
      obj.editCMap.Units = 'pixels';
      obj.editCMap.Position = [ 105 2 100 30];
      
      obj.ax.Units = 'pixels';
      xSize = max([5 (pixPos(3)-95)]);      
      ySize = max([5 0.95*(pixPos(4)-50)]);      
      obj.ax.Position = [70 50 xSize ySize];            
      
      %disp(['Setting TF Axes Position: ' num2str([70 50 xSize ySize])]);
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
        
    function updateImage(obj)
      axes(obj.ax); cla;
      set(crlEEG.gui.util.parentfigure.get(obj),'colormap',obj.cmap.cmap);
      obj.tfDecomp.imagesc('parent',obj.ax,...
                           'showBand',obj.showBand,...
                           'showChan',obj.showChan,...
                           'showTimes',obj.showTimes,...
                           'logImg',obj.logImg,...
                           'range',obj.imgRange,...
                           'colormap',obj.cmap)
                         
      if isempty(obj.cbar_)||~ishghandle(obj.cbar_)
        obj.cbar_ = colorbar('peer',obj.ax,'East');                   
      end;
      obj.cbar_.FontSize = 20;
      obj.cbar_.FontWeight = 'bold';
      
      tvals = obj.cbar_.Ticks;
      tickLabels = strsplit(num2str(obj.cmap.range(1) + tvals*(obj.cmap.range(2)-obj.cmap.range(1))));
      obj.cbar_.TickLabels = tickLabels;

      %obj.ax.Position = [0.025 0.05 0.965 0.9];
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
