classdef tfDecomp <  labelledArray
% Object class for time-frequency decompositions
%
%
% function obj = tfDecomp(decompType,tfData,tVals,fVals,chanLabels)
%
% Properties
% ----------
%   decompType : Name of the decomposition
%   tfData  : nFrequency X nTime x nChannel array of decomposition parameters
%   tVals   : Time values for each column in tfData
%   fVals   : Frequency values for each row in tfData
%  chanLabels: Cell string array of channel chanLabels
%
% Based on the labelledArray class available at:
%   https://github.com/hdamon/labelledArray
%
  
  properties (Dependent = true)
    decompType   
    dataType
    % Channel Labels
    chanLabels

    % Time Properties
    tfData
    tVals
    tRange

    % Frequency Properties
    fVals    
    fRange    
  end

  properties (Access=protected)
    decompType_
    dataType_
  end
  
  properties (Access=private,Constant)
    % Provided to avoid confusion elsewhere
    fDim = 1;
    tDim = 2;
    chanDim = 3;
  end

  properties
    decompParams % Used to store decomposition parameters
  end
  
  methods
    
    function obj = tfDecomp(tfData,varargin)
      %% Object Constructor Function
      
      % Configure standard dimensions
      dim(obj.fDim)    = arrayDim('dimName','frequency');
      dim(obj.tDim)    = arrayDim('dimName','time');
      dim(obj.chanDim) = arrayDim('dimName','channel');
      obj.dimensions = dim;
           
      if nargin>0
        p = inputParser;
        p.addParameter('decompType','',@ischar);
        p.addParameter('dataType','',@ischar);
        p.addParameter('tVals',[],@(x) isnumeric(x)&&isvector(x));
        p.addParameter('fVals',[],@(x) isnumeric(x)&&isvector(x));
        p.addParameter('chanLabels',[],@(x) ischar(x)||iscellstr(x));
        p.addParameter('decompParams',[]);
        p.parse(varargin{:});
        
        obj.tfData = tfData;
        obj.decompType = p.Results.decompType;
        obj.dataType   = p.Results.dataType;
        obj.tVals      = p.Results.tVals;
        obj.fVals      = p.Results.fVals;
        obj.chanLabels = p.Results.chanLabels;
        obj.decompParams = p.Results.decompParams;
               
      end;
    end
    

    function out = get.dataType(obj)
      out = obj.dataType_;
    end;
    
    function set.dataType(obj,val)
      obj.dataType_ = val;
    end;
    
    %% Get Time Range (No Set Method)
    function out = get.tRange(obj), out = [obj.tVals(1) obj.tVals(end)]; end;
    
    %% Get Frequency Range (No Set Method)
    function out = get.fRange(obj), out = [obj.fVals(1) obj.fVals(end)]; end;
    
    %% Set/Get Methods for obj.decompType
    % Just a redirect to the internal pro
    function out = get.decompType(obj), out = obj.decompType_; end;   
    function set.decompType(obj,val),   obj.decompType_ = val; end;
    
    %% Set/Get Methods for obj.chanLabels
    function out = get.chanLabels(obj)
      out = obj.dimLabels{obj.chanDim};
    end

    function set.chanLabels(obj,val)
      % Assign default channel values if 
      % provided with an empty array
      if isempty(val)
       nChan = obj.size(obj.chanDim);
       val = cell(nChan,1); 
       for i = 1:nChan
         val{i} =  ['Chan' num2str(i)];
       end
      end
      obj.dimLabels{obj.chanDim} = val; 
    end;
   
    %% Set/Get Methods for obj.tfData
    function out = get.tfData(obj)
      out = obj.array;
    end
    
    function set.tfData(obj,val)
      obj.array = val;
    end
   
    %% Set/Get Methods for obj.tVals
    function out = get.tVals(obj)
      out = obj.dimValues{obj.tDim};
    end
    
    function set.tVals(obj,val)
      obj.dimValues{obj.tDim} = val;
    end
    
    %% Set/Get Methods for obj.fVals
    function out = get.fVals(obj)
      % No default value for frequencies
      out = obj.dimValues{obj.fDim};
    end
   
    function set.fVals(obj,val)
      obj.dimValues{obj.fDim} = val;
    end
          
    function out = cat(dim,obj,a,varargin)      
    
      % Return 
      if ~exist('a','var'),a = [] ; end;
      if isempty(a), out = obj; return; end;
      if isempty(obj), out = a; return; end;
      
      assert(isEmptyOrEqual(obj.decompType,a.decompType),...
                  'Inconsistent decompType in concatenation');
      assert(isEmptyOrEqual(obj.dataType,a.dataType),...
                  'Inconsistent dataType in concatenation');
             
      out = cat@labelledArray(dim,obj,a);
      if isempty(out.decompType), out.decompType = a.decompType; end;
      if isempty(out.dataType), out.dataType = a.dataType; end;
              
      if ~isempty(varargin)
        if numel(varargin)<5
         out = cat(dim,out,varargin{:});
        else
         % Split recursion when concatenating multiple objects          
         nSplit = ceil(numel(varargin)/2);                
         blockA = cat(dim,out,varargin{1:nSplit});
         blockB = cat(dim,varargin{(nSplit+1):end});        
         out = cat(dim,blockA,blockB);
        end;        
        % Recurse when concatenating multiple objects
        %out = cat(dim,out,varargin{:});
      end;      
      
      function isValid = isEmptyOrEqual(A,B)
        isValid = isempty(A)||isempty(B);
        isValid = isValid||isequal(A,B);               
      end
      
    end;
    
    function out = subtract_baseline(obj,baseline)
      % Subtract a baseline frequency spectrum from all tfData columns.      
      assert(size(baseline,1)==size(obj.tfData,1),...
                'Incorrect Baseline Size');
      
      out = obj;
      out.tfData = abs(out.tfData) - repmat(baseline,1,size(out.tfData,2));
    end        
    
    function tfOut = selectTimes(tfIn,timesOut,varargin)
      % Select a subset of times from a time-frequency decomposition
      %
      % function tfOut = selectTimes(tfIn,timesOut)
      %
      % Inputs
      % ------
      %      tfIn : tfDecomp object
      %  timesOut : Timepoints to include in the output
      %               tfDecomp
      %
      % Outputs
      % -------
      %  tfOut : tfDecomp object with subselected times.
      %
      % Part of the crlEEG project
      % 2009-2018
      %
            
      % Output range must be sorted
      assert(issorted(timesOut),'Output times must be sorted');
      inRange = ( timesOut(1) >= tfIn.tVals(1) ) & ( timesOut(end) <= tfIn.tVals(end));
      
      % Just drop them?
      timesOut(timesOut<tfIn.tVals(1)) = [];
      timesOut(timesOut>tfIn.tVals(end)) = [];
      
      %assert(inRange,'Requested times are out of range');
      
      % Get indices of the time range to search in.
      [~,searchStart] = min(abs(tfIn.tVals-timesOut(1)));
      [~,searchEnd]   = min(abs(tfIn.tVals-timesOut(end)));
            
      if (searchEnd-searchStart)<numel(timesOut)
        outIdx = searchStart:searchEnd;
      else
      
      % Initialize output index
      outIdx = nan(numel(timesOut),1);
      outIdx(1) = searchStart;
      outIdx(end) = searchEnd;
                  
      idxSearch = searchStart;
      for idxOut = 2:(numel(timesOut)-1)
        minDeltaT = 1e100; % initialize minimum                
        deltaT = abs(tfIn.tVals(idxSearch)-timesOut(idxOut));
        while deltaT<minDeltaT          
          minDeltaT = deltaT;  % Found a new minimum
          idxSearch = idxSearch+1; % Advance to next time point
          if idxSearch > size(tfIn,2), break; end;
          deltaT = abs(tfIn.tVals(idxSearch)-timesOut(idxOut));  % Update deltaT       
        end;
        
        idxSearch = idxSearch-1; % Take one step back. 
        outIdx(idxOut) = idxSearch; % Update output                        
      end
      end;
      % Don't duplicate points in the output?
      outIdx = unique(outIdx);
      
      % Select the appropriate points to output.
      s.type = '()';
      s.subs = {':' outIdx ':'};
      
      tfOut = tfIn.subsref(s);
            
    end
    
    function p = show(obj,varargin)
      % Use the function from the .gui subpackage
      p = MatTSA.gui.tfDecomp.showTF(obj,varargin{:});
    end
    
    function varargout = imagesc(obj,varargin)
      % Overloaded imagesc method for tfDecomp objects
      %
      % Inputs
      % ------
      %
      % Optional Inputs
      % ------
      %  range: Range to display
      %
      % Optional Param-value Inputs
      % ---------------------------
      %  showChan : Index of channel to display
      %    logImg : Flag to turn on logarithmic display
      %  showBand : 1x2 Array: [lower upper] frequency band
      %    parent : Matlab gui handle to parent to 
      %
      %
                      
      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('range',[],...
              @(x) isempty(x)||(isnumeric(x)&&isvector(x)&&(numel(x)==2)));
      p.addParameter('showChan',1);
      p.addParameter('logImg',false);
      p.addParameter('showBand' , [] , @(x) isempty(x)||(isnumeric(x)&&isvector(x))); 
      p.addParameter('showTimes', [] , @(x) isempty(x)||(isnumeric(x)&&isvector(x))); 
      p.addParameter('parent'   , [] , @(x) ishghandle(x));
      p.addParameter('colormap',guiTools.widget.alphacolor,@(x) isa(x,'guiTools.widget.alphacolor'));
      p.parse(varargin{:});
      
      %% Get Channels to Show
      if isempty(p.Results.showChan)
        idxChan = 1;
      else
        idxChan = p.Results.showChan;
      end;
      

      %% Get Frequency Band Indices
      if ~isempty(p.Results.showBand)
        [~,idxLow] = min(abs(obj.fVals-p.Results.showBand(1)));
        [~,idxHi ] = min(abs(obj.fVals-p.Results.showBand(2)));
        idxF = idxLow:idxHi;
      else
        idxF = 1:numel(obj.fVals);
      end
      
      %% Get Time Indices
      if ~isempty(p.Results.showTimes)
        if numel(p.Results.showTimes)==2
          % Treat it as a range
          [~,idxLow] = min(abs(obj.tVals-p.Results.showTimes(1)));
          [~,idxHi ] = min(abs(obj.tVals-p.Results.showTimes(2)));
          idxT = idxLow:idxHi;
        else
          idxT = p.Results.showTimes;
        end;
      else
        idxT = ':';
      end;
      
      %% Get the Image to Display.
      s(1).type = '.';
      s(1).subs = 'tfData';
      s(2).type = '()';
      s(2).subs = {idxF idxT idxChan};
      
      showImg = obj.subsref(s);
      s(2).subs = {idxF ':' idxChan};
      rangeImg = obj.subsref(s);
           
      % Take Magnitude if Complex.
      if any(any(imag(showImg)))
        showImg = abs(showImg);
        rangeImg = abs(rangeImg);
      end;
                    
      %% Get Image Range
      if isempty(p.Results.range)
        imgRange(1) = prctile(rangeImg(:),0.001);
        imgRange(2) = prctile(rangeImg(:),99.999);
      else
        imgRange = p.Results.range;
      end

      %% If log requested
      if p.Results.logImg
        showImg = log10(showImg);        
        imgRange = log10(imgRange);                
      end;                  
      
      %% Get the RGB Image
      cmap = p.Results.colormap;    
      if isempty(cmap.range)||isequal(cmap.range,[0 1])
        % Only override if it's the default
        cmap.range = imgRange;
      end;      
                  
      [rgb,alpha] = cmap.img2rgb(showImg);
      tData = obj.tVals(idxT);
      fData = obj.fVals(idxF);

      %% Get Object Parent
      if isempty(p.Results.parent)
        par = figure;
      else
        if ishghandle(p.Results.parent,'figure')
         figure(p.Results.parent);
        elseif ishghandle(p.Results.parent,'axes')
          axes(p.Results.parent);
        end
      end;
       
      yyaxis left;
      set(gca,'YColor',[0 0 0]);
      img = image(tData,[],rgb,'AlphaData',alpha);
      
      currAxis = gca;
      showF = obj.fVals(idxF(floor(currAxis.YTick)));
      for i = 1:numel(showF)
        fLabel{i} = num2str(showF(i));
      end;
      currAxis.YTickLabel = fLabel;
      
      
      set(gca,'YDir','normal');
      ylabel('Frequency');
      xlabel('Time');
      
      if nargout>0
        varargout{1} = img;
      end;
      
    end % imagesc()
      
    function out = PSD(obj,varargin)
      % Convert a time-frequency decomposition to power spectral density                  
      out = obj.copy;
      if numel(varargin)>0
        s.type = '()';
        s.subs = {varargin{:}};
        out = subsref(out,s);
      end
      if ~isequal(obj.dataType,'PSD')      
        out.tfData = abs(out.tfData).^2;
        out.decompType = [obj.decompType '_PSD'];
        out.dataType = 'PSD';      
      end;
    end
    
    function out = phase(obj,varargin)
      % Returns the phase of a complex valued tfDecomp.
      %
      
      out = obj.copy;
      if numel(varargin)>0
        s.type = '()';
        s.subs = {varargin{:}};
        out = subsref(out,s);
      end;
      
      X = real(out.data);
      Y = imag(out.data);
      
      phase = atan(Y./abs(X));
      phase((X==0)&(Y==0)) = nan;
            
      out.tfData = phase;
      out.decompType = [obj.decompType '_Phase'];
      
      
    end
    
    function out = abs(obj)
      % Convert a time-frequency decomposition to spectral magnitude            
      out = obj.copy;
      out.tfData = abs(obj.tfData);
      out.decompType = [obj.decompType '_ABS'];      
    end;
    
    function out = PLF(obj)  
      % Convert a time-frequency decomposition to Phase Locking Factor
      %
      % (NEEDS TO BE AVERAGED ACROSS A LOT OF DECOMPOSITIONS)            
      out = obj.copy;
      out.decompType = [obj.decompType '_PLF'];
      tmp = obj.tfData;
      out.tfData = tmp./abs(tmp);
      out.dataType = 'PLF';
    end;
              
    tfOut = interpFrequencies(tfIn,outFreq,type);
    
    
  end
  
  %% Protected Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Access=protected)
    
        %% SubCopy
    function out = subcopy(obj,varargin)
      % Copy object, including only a subset of timepoints and columns. If
      % not provided or empty, indices default to all values.
      %
     
      [out, chanIdx] = obj.subcopy@labelledArray(varargin{:});                          
    end
    
  end
  
end

