classdef timeseries < labelledArray
  % Data class for timeseries data  
  %
  % obj = MatTSA.timeseries(data,chanLabels,varargin)
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
    data % Redirected from obj.a
  
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
        
  properties (Access=private, Constant)
    % Useful to prevent confusion elsewhere
    tDim    = 1;
    chanDim = 2;
  end
  
  methods
    
    function obj = timeseries(varargin)
      
      %% Input Parsing
      obj = obj@labelledArray;
      
      dims(obj.tDim) = arrayDim('dimName','time');
      dims(obj.chanDim) = arrayDim('dimName','channel');      
      obj.dimensions = dims;
      
      if nargin>0
        p = inputParser;
        p.addRequired('data',@(x) (isnumeric(x)&&ismatrix(x))||...
                                    isa(x,'MatTSA.timeseries'));
        p.addParameter('chanLabels', []  ,@(x) isempty(x)||iscellstr(x));
        p.addParameter(  'chanType', []  ,@(x) ischar(x)||iscellstr(x));
        p.addParameter(     'tVals', []  ,@(x) isempty(x)||isvector(x));
        p.addParameter(    'tUnits','sec',@ischar);
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
          units = repmat(units,numel(label),1);
        end;
        
        if ~iscellstr(type)
          type = repmat(typs,numel(label),1);
        end;
        
        for i = 1:numel(label)
          addChannel(obj,data(:,i),label{i},units{i},type{i});
        end;
        
        return;
      end
      
      % Defaults
      if ~exist('units','var')||isempty(units), units = '_'; end;
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
         obj.array_ = [obj.data data(:)];        
         obj.dimLabels{obj.chanDim}{end+1} = label;
         obj.dimUnits_{obj.chanDim}{end+1} = units;
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
      obj.dimLabels_{obj.chanDim} = obj.chanLabels(idx);
      obj.dimUnits_{obj.chanDim} = obj.tUnits(idx);
      obj.array_ = obj.data(:,idx);
      obj.chanType_ = obj.chanType_(idx);
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
      obj.updateDataRange;
    end % END set.chanType
                
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
      obj.updateDataRange;
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
      out = obj.dimValues{obj.tDim};      
    end    
    function set.tVals(obj,val)
      if isempty(val)
        % Default time values.
        val = (1./obj.sampleRate)*(0:size(obj.data,1)-1);
      end
      obj.dimValues{obj.tDim} = val;
    end;
    
    %% Get/Set Methods for obj.tUnits
    function out = get.tUnits(obj)
      out = obj.dimUnits{obj.tDim};
    end
    function set.tUnits(obj,val)
      obj.dimUnits{obj.tDim} = val;
    end;
    
    
    %% Get/Set Methods for obj.dataRange
    function rangeOut = get.dataRange(obj)
      rangeOut = obj.dataRange_;
%       dataChans = obj.getChannelsByType('data');   
%       if ~any(dataChans)
%         rangeOut = [0 1]; return;
%       end;
%       rangeOut = [min(min(obj.data(:,dataChans))) ...
%                   max(max(obj.data(:,dataChans)))];                
    end;                  
    function set.dataRange(obj,~)
      error('obj.dataRange is derived from obj.data');
    end;
        
    function updateDataRange(obj)
      dataChans = obj.getChannelsByType('data');   
      if ~any(dataChans)
        obj.dataRange_ = [0 1]; return;
      end;
      obj.dataRange_ = [min(min(obj.data(:,dataChans))) ...
                       max(max(obj.data(:,dataChans)))];  
    end
    
    %% Get/Set Methods for obj.tRange
    function rangeOut = get.tRange(obj)
      rangeOut = [obj.tVals(1) obj.tVals(end)];
    end;    
    function set.tRange(obj,~)
      error('obj.tRange is derived from obj.tVals');
    end;     
    
    %% Methods with their own m-files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotOut = butterfly(tseries,varargin);
    plotOut = plotWithDecomp(tseries,varargin);
    
    function outTseries = filtfilt(tseries,dFilter)
      % Overloaded filtfilt function for crltseries.type.timeseries objects
      %
      % Inputs
      % ------
      %   tseries : A crltseries.type.tseries.object to be filtered
      %  dFilter : A Matlab digital filter (typically created with designfilt)
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
    function out = subcopy(obj,varargin)
      % Not really necessary, but included for the future.                
      [out,dimIdx] = obj.subcopy@labelledArray(varargin{:});       
      out.chanType_ = out.chanType_(dimIdx{2});
      
      if ~isempty(obj.decomposition)
        decompNames = fields(obj.decomposition);
        for i = 1:numel(decompNames)
          tmpDecomp = obj.decomposition.(decompNames{i}).copy;
          tmpDecomp = tmpDecomp.selectTimes(out.xvals);
          out.decomposition.(decompNames{i}) = tmpDecomp;
        end
        
      end      
    end
  end
  
  %% Static Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Static=true)
  end
  
end
