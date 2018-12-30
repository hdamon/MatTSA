classdef timeseries < labelledArray
  % Data class for timeseries data  
  %
  % obj = MatTSA.timeseries(data,chanLabels,varargin)
  %
  % Subclass of labelledArray that adds properties and methods suitable for
  % general time series analysis.
  %
  % Properties:
  % -----------------
  %   data:
  %  chanLabels:  
  %     chanType : 
  %  tUnits
  %  tVals
  %  tRange
  %
  %   
  %
  % Inputs
  % ------
  %   data : nTime x nChannels array of time series data
  %   chanLabels : (Optional) Cell array of length nChannels containing label strings
  %   
  % Param-Value Pairs
  % -----------------
  %     tUnits : Units for the data (DEFAULT: [])
  %  dataUnits : Units of time (DEFAULT: 'sec')
  %      tVals : Timings associated with each sample. 
  % sampleRate : Sample rate for the data (DEFAULT: 1Hz)
  % 
  % Referencing into timeseries objects
  % -----------------------------------
  %   As a subclass of labelledArray, MatTSA.timeseries has the same
  %   referencing options available. "help matTSA.timeseries" for more
  %   information.
  % 
  %
  % Written By: Damon Hyde
  % Part of the MatTSA Package
  % 2018-
  %
  
  
  properties (Dependent = true)
    data % Redirected from obj.array
  
    % Channel Parameters
    chanLabels
    chanType

    % Time Parameters:
    tUnits
    tVals  
    tRange    
    sampleRate

    % Data Parameters
    dataUnits
    dataRange
  end;
    
  properties
    decomposition
  end
  
  properties (Access=protected)   
    sampleRate_;   
    chanType_;
    chanInfo_;
    dataRange_;
  end;
        
  properties (Access=protected, Constant)
    % Useful to prevent confusion elsewhere
    timeDim = 1;
    chanDim = 2;
  end
  
  methods
    
    function obj = timeseries(varargin)
      
      %% Input Parsing
      obj = obj@labelledArray;
      
      dims(obj.timeDim) = arrayDim('dimName','time');
      dims(obj.chanDim) = arrayDim('dimName','channel');      
      obj.dimensions = dims;
      
      if nargin>0
        p = inputParser;
        p.addRequired('data',@(x) (isnumeric(x)&&(ndims(x)>=2))||...
                                    isa(x,'MatTSA.timeseries'));
        p.addParameter('chanLabels', []  ,@(x) isempty(x)||iscellstr(x));
        p.addParameter(  'chanType', []  ,@(x) ischar(x)||iscellstr(x));
        p.addParameter(     'tVals', []  ,@(x) isempty(x)||isvector(x));
        p.addParameter(    'tUnits',[],@ischar);
        p.addParameter('sampleRate',  1  ,@(x) isnumeric(x)&&isscalar(x));
        p.addParameter( 'dataUnits', []  ,@(x) ischar(x)||iscellstr(x));
                        
        p.parse(varargin{:});
        
        if isa(p.Results.data,'MatTSA.timeseries')
          % Return a new object with the same values as the input
          obj = obj.copyValuesFrom(p.Results.data);
          return;
        end                
        
        %% Set Object Properties
        obj.data       = p.Results.data;
        obj.chanLabels = p.Results.chanLabels;
        obj.chanType   = p.Results.chanType;
        obj.tVals      = p.Results.tVals;
        obj.tUnits     = p.Results.tUnits;
        obj.sampleRate = p.Results.sampleRate;
        obj.dataUnits  = p.Results.dataUnits;
        
      end;
    end         
    
    output = blockProcess(tSeriesIn,procHandle,varargin);
    
    %% Main MatTSA.timeseries plotting function
    function out = plot(obj,varargin)
      % Overloaded plot function for MatTSA.timeseries objects
      %
      % Inputs
      % ------
      %   obj : MatTSA.timeseries object
      % 
      % Param-Value Pairs
      % -----------------
      %   'type' : Type of plot to display (DEFAULT: 'dualplot')
      %              Valid Options:
      %                'dualplot'  : 
      %                'butterfly' : Plot each channel on top of one
      %                               another in the same axis.
      %                'split'     : Display each channel as it's own
      %                               plot in the same axis
      %               
      %
      %
      % All other inputs are passed directly to the plotting function.
      %
      
      p = inputParser;
      p.KeepUnmatched = true;
      p.addParameter('type','dualplot',@(x) ischar(x));
      p.parse(varargin{:});
            
      switch lower(p.Results.type)
        case 'dualplot'
          out = MatTSA.gui.timeseries.dualPlot(obj,p.Unmatched);
        case 'withdecomp'
          out = plotWithDecomp(obj,p.Unmatched);
        case 'butterfly'
          out = butterfly(obj,p.Unmatched);
        case 'split'
          out = split(obj,p.Unmatched);
        otherwise
          error('Unknown plot type');
      end                     
    end

    %% Plotdata is modified a bit to improve display
    function out = getPlotData(obj)      
      % Return modified timeseries data for plotting
      %
      % out = getPlotData(obj)
      %
      % Data channels are returned unchanged. Boolean and auxillary
      % channels have their amplitude modified so that they are in the
      % range [0 max(obj.data(:,<datachannels>)]
      %
      % This is intended to help with plot scaling when 
      % plotting multiple channels.
      %
      out = obj.data;
      
      % Boolean Data Channels
      boolChan = obj.getChannelsByType('bool');
      if ~isempty(boolChan)
        out(:,boolChan) = 0.75*obj.dataRange(2)*out(:,boolChan); 
      end;
      
      % Auxilliary Channels
      auxChan = obj.getChannelsByType('aux');      
      for i = 1:numel(auxChan)        
          m = max(abs(out(:,auxChan(i))));
          out(:,auxChan(i)) = 0.75*(obj.dataRange(2)/m)*out(:,auxChan(i));        
      end
        
    end    
    
    %% Add/Remove Channels from a Timeseries Objects
    function addChannel(obj,data,label,units,type,replace)
      % Add a channel to a timeseries object
      %
      % function addChannel(obj,data,label,units,type,replace)
      %
      % Inputs
      % ------
      %   obj   : Timeseries object to add channel to
      %  data   : Vector containing 
      %  label  : Channel label(s)
      %  units  : Physical Units
      %   type  : Type of Channel ('data','aux','bool')
      % replace : When set to true, replaces an existing channel
      %
      % If the input data is a matrix, label must be a cell string of
      % channel chanLabels. Units and type can then either be a single
      % character string (Uniform across channels), or cell arrays with
      % individual values for each channel.
      % 
        
      if ~exist('replace','var'), replace = false; end;
      if replace
        obj.removeChannel(label);
      end
              
      if (size(data,2)==size(obj,1))&&(size(data,1)~=size(obj,1))
        data = data';
      end;
      
      assert(size(data,1)==size(obj,1),'Channel Data Size is Incorrect');
      test1 = iscellstr(label)&&(size(data,2)==numel(label));
      test2 = size(data,2)==1;
      assert(test1||test2,'Incorrect number of chanLabels provided');
                  
      % Recurse
      if iscellstr(label)
        if ~exist('units','var')||isempty(units), units = repmat({'_'},numel(label),1); end;
        if ~exist('type','var')||isempty(type), type = repmat({'data'},numel(label),1); end;
        
        if ~iscellstr(units)
          [tmp{1:numel(label)}] = deal(units);
          units = tmp;          
        end;
        
        if ~iscellstr(type)
          [tmp{1:numel(label)}] = deal(type);
          type = tmp;          
        end;
        
        for i = 1:numel(label)
          addChannel(obj,data(:,i),label{i},units{i},type{i});
        end;
        
        return;
      end
      
      % Defaults
      if ~exist('units','var')||isempty(units), units = []; end;
      if ~exist('type','var')||isempty(type), type = 'data'; end;
      if ~exist('replace','var'), replace = false; end;
      
      % Add a single label
      if ismember(label,obj.chanLabels)
        if replace
          warning('Channel replacement unimplemented');
          return;
        else
          error('Channel already exists');
        end;
      else
         %tmpTSeries = MatTSA.timeseries(data,'chanLabels',label,'chanUnits',units
        
        
         obj.array_ = cat(2,obj.data,data);        
         
         tmpDim = arrayDim;
         tmpDim.dimName = obj.dimensions(obj.chanDim).dimName;
         tmpDim.dimLabels = {label};         
         if isempty(units)
           tmpDim.dimUnits = [];
         else
           tmpDim.dimUnits = {units};
         end;
         newDim = cat(1,obj.dimensions(obj.chanDim),tmpDim);
         obj.dimensions(obj.chanDim) = newDim;                           
         obj.chanType_{end+1} = type;         
      end
    end;
        
    function removeChannel(obj,label)
      % Remove one or more channels from a timeseries object
      %
      % function removeChannel(obj,label)
      %
      % Inputs
      %    obj : MatTSA.timeseries object
      %  label : List of channel chanLabels to remove. Can be either a string,
      %           or a cell array of strings.
      %
      
      if ~iscell(label), label = {label}; end;            
      assert(iscellstr(label),'Labels must be provided as strings');      
      idx = ~ismember(obj.chanLabels,label);
           
      % Truncate the internal channels
      obj.array_ = obj.data(:,idx);
      obj.dimensions(obj.chanDim) = ...
                      obj.dimensions(obj.chanDim).subselectDimensions(idx);                  
      obj.chanType_ = obj.chanType_(idx);
    end
           
    function out = cat(dim,obj,a,varargin)
      % Concatenate timeseries objects
      %
      %
      
      assert(isa(a,class(obj)),'Can only concatenate like objects');
     
      
      
      out = cat@labelledArray(dim,obj,a);
                  
      % Pick correct sample rate
        if isMostlyEmpty(obj.dimensions(1))
          out.sampleRate = a.sampleRate;
        elseif isMostlyEmpty(a.dimensions(1))
          out.sampleRate = obj.sampleRate;
        else
          assert(a.sampleRate==obj.sampleRate,'Mismatched sample rates');
        end;      
      
      if ~isempty(varargin)
        % Recurse when concatenating multiple objects
        out = cat(dim,out,varargin{:});
      end;
      
    end    
    
    %% Retrieve Channels By Type
    function out = isChannelType(obj,val)
      % Returns a logical array that is true if a channels dataUnits type
      % matches val
      if isempty(obj.chanType)
        out = false(size(obj,2),1);
      else
        out = cellfun(@(x) isequal(x,val),obj.chanType);
      end
    end
    
    function out = getChannelsByType(obj,val)
      out = find(obj.isChannelType(val));
    end
    
    %% Retrieve Channels By Physical Units
    function out = isUnitType(obj,val)
      % 
      % function out = isUnitType(obj,val)
      %
      % Returns a logical array 
      out = cellfun(@(x) isequal(x,val),obj.dataUnits);
    end
    
    function out = getChannelBdataUnits(obj,val)
      out = find(obj.isUnitType(obj,val));
    end
            
    %% GET/SET METHODS FOR DEPENDENT PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Get/Set Methods for obj.chanType
    function out = get.chanType(obj)      
      if isempty(obj.chanType_)
        % Initialize to defaults if empty
        obj.chanType = [];
      end;
      out = obj.chanType_;      
    end; % END get.chanType    
    
    function set.chanType(obj,val)
      if isempty(val), val = 'data'; end;
      assert(ischar(val)||iscellstr(val),...
              'chanType must be a character string or cell array of strings');
      if ~iscellstr(val)
        [cellVal{1:size(obj,2)}] = deal(val); 
      else
        cellVal = val;
      end;
      
      assert(numel(cellVal)==size(obj,2),...
              'chanType must have a number of elements equal to the number of channels');
      obj.chanType_ = cellVal;
      obj.arrayRange_ = [];      
    end % END set.chanType
                
    function set.chanType_(obj,val)
      if size(val,1)>1
        val = val';
      end;
      obj.chanType_ = val;
    end
    
    %% Get/Set Methods for obj.dataUnits
    function out = get.dataUnits(obj)
      out = obj.dimUnits{obj.chanDim};      
    end;    
    
    function set.dataUnits(obj,val)
      if isempty(val), obj.dimUnits{obj.chanDim} = []; return; end;
      assert(ischar(val)||iscellstr(val),...
              'tUnits must be a character string or cell array of strings');
      if ~iscellstr(val)
        [cellVal{1:size(obj,2)}] = deal(val); 
      else
        cellVal = val;
      end;
            
      obj.dimUnits{obj.chanDim} = cellVal;                        
    end    
    
    %% Get/Set Methods for Data
    function out = get.data(obj)
      out = obj.array;
    end    
    
    function set.data(obj,val)
      obj.array = val;      
    end
           
    %% Set/Get Methods for obj.chanLabels
    function out = get.chanLabels(obj)           
      out = obj.dimLabels{obj.chanDim};      
    end % END get.chanLabels
    
    function set.chanLabels(obj,val)
      if isempty(val) 
        val = cell(1,size(obj.data,2));
        for i = 1:size(obj.data,2)
          val{i} = ['Chan' num2str(i)];
        end      
      end;
      
      obj.dimLabels{2} = val;      
    end % END set.chanLabels
            
    %% Get/Set Methods for obj.sampleRate
    function out = get.sampleRate(obj)
      if ~isempty(obj.sampleRate_)
       out = obj.sampleRate_;
      else
       out = 1;
      end;
    end;   
    function set.sampleRate(obj,val)
      if isempty(val), obj.sampleRate = []; return; end;
      assert(isnumeric(val)&&isscalar(val),...
         'Sample rate must be a scalar numeric value');
       obj.sampleRate_ = val;
    end
            
    %% Get/Set Methods for obj.tVals    
    function out = get.tVals(obj)
      out = obj.dimValues{obj.timeDim};      
    end    
    function set.tVals(obj,val)
      if isempty(val)
        % Default time values.
        val = (1./obj.sampleRate)*(0:size(obj.data,1)-1);
      end
      obj.dimValues{obj.timeDim} = val;
      obj.postTValUpdate;
    end;   
    
    function postTValUpdate(obj)
      % Function called after obj.tVals are updated
      return        
    end
    
    %% Get/Set Methods for obj.tUnits
    function out = get.tUnits(obj)
      out = obj.dimUnits{obj.timeDim};
    end
    function set.tUnits(obj,val)
      obj.dimUnits{obj.timeDim} = val;
    end;
    
    
    %% Get/Set Methods for obj.dataRange
    function rangeOut = get.dataRange(obj)
      %% Should be updated to use the arrayRange functionality of labelledArray
      rangeOut = obj.arrayRange;
    end;                  
        
    function updateArrayRange(obj,force)
      % Timeseries objects only consider 'data' channels in computing the
      % arrayRagne
      if ~exist('force','var'),force = false; end;
      if isempty(obj.dataRange_)&&~force      
        return;
      end;
      
      dataChans = obj.getChannelsByType('data');
      if ~any(dataChans)
        obj.dataRange_ = [0 1]; return;
      end;
      obj.updateArrayRange@labelledArray(force,{':',dataChans});
    end
    
    %% Get/Set Methods for obj.tRange
    function rangeOut = get.tRange(obj)
      rangeOut = obj.dimensions(obj.timeDim).dimRange;
    end;    
    function set.tRange(obj,~)
      error('obj.tRange is derived from obj.tVals');
    end;     
    
    function [out,varargout] = applyDimFunc(funcHandle,obj,idxDim,varargin)
      if nargout==1
        out = applyDimFunc@labelledArray(funcHandle,obj,idxDim,varargin{:});
      else
        outCell = cell(1,nargout-1);
        [out,outCell{:}] = applyDimFunc@labelledArray(funcHandle,obj,idxDim,varargin{:});
        varargout = outCell;
      end
      if varargin{idxDim}==2
        % New Channel Type if Collapsing Channel Dimension
       for i = 1:numel(out)
         out(i).chanType = {'data'};
       end;
      end;
    end
    
    %% Methods with their own m-files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotOut = butterfly(tseries,varargin);
    plotOut = plotWithDecomp(tseries,varargin);
    
    function outTseries = filtfilt(tseries,dFilter)
      % Overloaded filtfilt function for MatTSA.timeseries objects
      %
      % outTseries = filtfilt(tseries,dFilter)
      %
      % Inputs
      % ------
      %   tseries : A crltseries.type.tseries.object to be filtered
      %  dFilter : A Matlab digital filter (typically created with designfilt)
      %
      % Output
      % ------
      %   outTseries : A new MatTSA.timeseries object, copied from the
      %                 original, with the specified filter applied to
      %                 all data channels.
      %
                  
      dataChans = tseries.getChannelsByType('data');      
      tmp = filtfilt(dFilter,tseries.data(:,dataChans));      
      outTseries = tseries.copy;
      outTseries.data(:,dataChans) = tmp;      
    end
    
    function outTseries = filter(tseries,dFilter)
      % Overloaded filter function for MatTSA.timeseries objects
      %
      % Inputs
      % ------
      %  tseries : A MatTSA.timeseries object to be filtered
      % dFilter : A matlab digitalFilter object (typically created with
      % designfilt)
      %      
      tmp = filter(dFilter,tseries.data(:,tseries.getChannelsByType('data')));
      
      outTseries = tseries.copy;
      outTseries.data(:,outTseries.getChannelsByType('data')) = tmp;      
    end
        
  end
  
  methods (Access=protected)
     function obj = copyValuesFrom(obj,valObj)
      % Individually copy values from a second object
      obj = obj.copyValuesFrom@labelledArray(valObj);      
      if isa(valObj,'MatTSA.timeseries')
       % Can only copy these if it's actually a timeseries object.
       obj.sampleRate = valObj.sampleRate;    
      end;
    end
    
    %% SubCopy
    function [out,varargout] = subcopy(obj,varargin)
            
      % Subselection the Array and Dimensions
      [out,dimIdx] = obj.subcopy@labelledArray(varargin{:});       
      
      % Subselect Channel Types
      out.chanType_ = obj.chanType(dimIdx{2});
      
      % Subselect any Decompositions
      if ~isempty(obj.decomposition)
        decompNames = fields(obj.decomposition);
        for i = 1:numel(decompNames)
          tmpDecomp = obj.decomposition.(decompNames{i}).copy;
          tmpDecomp = tmpDecomp.selectTimes(out.tVals);
          out.decomposition.(decompNames{i}) = tmpDecomp;
        end        
      end  
      
      if nargout>1
        varargout{1} = dimIdx;
      end;
      
    end
  end
  
  %% Static Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Static=true)
  end
  
end
