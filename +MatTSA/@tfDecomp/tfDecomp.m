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
    
    function obj = tfDecomp(decompType,tfData,tVals,fVals,chanLabels)
      %% Object Constructor Function
      if nargin>0
        obj.decompType_ = decompType;
        obj.tfData = tfData;
        obj.tVals  = tVals;
        obj.fVals  = fVals;      
        obj.chanLabels = chanLabels;
      end;
    end
    

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
      obj.dimLabels{obj.chanDim} = l; 
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
      out = obj.dimValues{tDim};
    end
    
    function set.tVals(obj,val)
      obj.dimValues{tDim} = val;
    end
    
    %% Set/Get Methods for obj.fVals
    function out = get.fVals(obj)
      % No default value for frequencies
      out = obj.dimValues{fDim};
    end
   
    function set.fVals(obj,val)
      obj.dimValues{fDim} = val;
    end
          
              
    %% SubCopy
    function out = subcopy(obj,varargin)
      % Copy object, including only a subset of timepoints and columns. If
      % not provided or empty, indices default to all values.
      %
      % Mostly intended as a utility function to simplify subsref.
      %
      if ~exist('fIdx','var'), fIdx = ':'; end;
      if ~exist('tIdx','var'), tIdx = ':'; end;
      if ~exist('chanIdx','var'),chanIdx = ':'; end;      
      

        
      out = obj.copy;                 
      out.tVals_ = out.tVals_(tIdx);
      out.fVals_ = out.fVals_(fIdx);
      out.chanLabels_ = out.chanLabels(chanIdx);
      out.tfData_  = out.tfData_(fIdx,tIdx,chanIdx);        
      
    end
    
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
          deltaT = abs(tfIn.tVals(idxSearch)-timesOut(idxOut));  % Update deltaT       
        end;
        
        idxSearch = idxSearch-1; % Take one step back. 
        outIdx(idxOut) = idxSearch; % Update output                        
      end
      
      % Don't duplicate points in the output?
      %outIdx = unique(outIdx);
      
      % Select the appropriate points to output.
      s.decompType = '()';
      s.subs = {':' outIdx ':'};
      
      tfOut = tfIn.subsref(s);
            
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
      
      import crlEEG.util.validation.*
                  
      p = inputParser;
      p.KeepUnmatched = true;
      p.addOptional('range',[],@(x) emptyOk(x,@(x) isNumericVector(x,2)));
      p.addParameter('showChan',1);
      p.addParameter('logImg',false);
      p.addParameter('showBand',[],@(x) emptyOk(x,@(x) isNumericVector(x))); 
      p.addParameter('showTimes',[],@(x) emptyOk(x,@(x) isNumericVector(x))); 
      p.addParameter('parent',[],@(x) ishghandle(x));
      p.addParameter('colormap',crlEEG.gui.widget.alphacolor,@(x) isa(x,'crlEEG.gui.widget.alphacolor'));
      p.parse(varargin{:});
      
      %% Get Channels to Show
      if isempty(p.Results.showChan)
        idxChan = 1;
      else
        idxChan = p.Results.showChan;
      end;
      
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
      s(1).decompType = '.';
      s(1).subs = 'tfData';
      s(2).decompType = '()';
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
      
            
      img = image(tData,[],rgb,'AlphaData',alpha);
      
      currAxis = gca;
      showF = obj.fVals(idxF(round(currAxis.YTick)));
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
    
    function out = subsrefTFX(obj,varargin)
      % Not sure this is entirely needed. Might be a bit of a hack
      if ~isempty(varargin)
        s(1).decompType = '.';
        s(1).subs = 'tfData';
        s(2).decompType = '()';
        s(2).subs = varargin;
        out = obj.subsref(s);
      else
        out = obj.tfData;
      end             
    end
    
    function out = PSD(obj,varargin)
      % Convert a time-frequency decomposition to power spectral density                  
      out = obj.copy;
      out.tfData = abs(obj.subsrefTFX(varargin{:})).^2;
      out.decompType = [obj.decompType '_PSD'];
    end
    
    function out = abs(obj,varargin)
      % Convert a time-frequency decomposition to spectral magnitude            
      out = obj.copy;
      out.tfData = abs(obj.subsrefTFX(varargin{:}));
      out.decompType = [obj.decompType '_ABS'];      
    end;
    
    function out = PLF(obj,varargin)  
      % Convert a time-frequency decomposition to Phase Locking Factor
      %
      % (NEEDS TO BE AVERAGED ACROSS A LOT OF DECOMPOSITIONS)      
      out = obj.copy;
      out.decompType = [obj.decompType '_PLF'];
      tmp = obj.subsrefTFX(varargin{:});
      out.tfData = tmp./abs(tmp);
    end;
              

    
    function out = sqrt(obj)
      out = obj.copy;
      out.tfData = sqrt(out.tfData);
      out.decompType = ['sqrt(' out.decompType ')'];
    end
    
    function isValid = isConsistent(obj,b)
      % Check consistency between tfDecomp objects
      %
      % Used to check if math operations can be applied between two
      % timefrequency decompositions
      %
      isValid = isa(b,'tfDecomp');      
      if ~isValid, return; end;
            
      fValsEqual = true;
      if size(obj,1)==size(b,1)
       % If frequency dimensions are equal, check that they're the same
       fValsEqual = isequal(obj.fVals,b.fVals);
      end;
      
      tValsEqual = true;
      if size(obj,2)==size(b,2)
        % If time dimensions are equal, check that they're the same.
        tValsEqual = isequal(round(obj.tVals,8),round(b.tVals,8));
      end
      
      chanEqual = true;
%       if size(obj,3)==size(b,3)
%         chanEqual = isequal(obj.chanLabels,b.chanLabels);
%       end
      
      sizeValid = crlEEG.util.validation.arraySizeForBSXFUN(size(obj),size(b));
      
      isValid = sizeValid && tValsEqual && fValsEqual && chanEqual;      
    end
    
    function [fVals,tVals,chan] = consistentDimensions(obj,b)
      assert(isa(b,'tfDecomp'),...
                'Second input must be a tfDecomp object');
              
      assert(isConsistent(obj,b),'Inconsistent decomposition sizes');  
      
      sizeEqual = crlEEG.util.validation.compareSizes(size(obj),size(b));
      
      fields = {'fVals' 'tVals' 'chanLabels'};        
      fVals = [];
      tVals = [];
      chan = [];
      for i = 1:numel(sizeEqual)
        if sizeEqual(i)
          tmp = obj.(fields{i});
        else
          if size(obj,i)==1
            tmp = b.(fields{i});
          elseif size(b,i)==1
            tmp = obj.(fields{i});
          else
            error('Shouldn''t be getting here');
          end;
        end
        switch i
          case 1
            fVals = tmp;
          case 2
            tVals = tmp;
          case 3
            chan = tmp;
        end
      end                            
    end
    
    function decompOut = applyFcnToTFX(obj,b,funcHandle)
      % Use bsxfun to apply funcHandle to obj.tfData
      %
      % decompOut = applyFcnToTFX(obj,b,funcHandle)
      %
      % Inputs
      % ------
      %        obj : tfDecomp object
      %          b : Value to apply
      % funcHandle : Function handle
      %
      % Checks the consistency of the inputs, and then calls:
      %
      % tfData = bsxfun(funcHandle,obj,tfData,coeff)
      %
      % When b is:
      %   A tfDecomp obj:  coeff = b.tfData
      %                          OTHERWISE:  coeff = b;
      %
      %      
      
      for idxObj = 1:numel(obj)
      
      if isa(b,'tfDecomp')
        assert(isConsistent(obj(idxObj),b),'Inconsistent decomposition sizes');
        coeff = b.tfData;
        [fVals,tVals,chan] = consistentDimensions(obj(idxObj),b);
      else
        coeff = b;
        fVals = obj(idxObj).fVals;
        tVals = obj(idxObj).tVals;
        chan = obj(idxObj).chanLabels;
      end;
            
      tfData = bsxfun(funcHandle,obj.tfData,coeff);
      newType = [obj.decompType '_' func2str(funcHandle)];
      decompOut(idxObj) = tfDecomp(newType,tfData,tVals,fVals,chan);      
      
      end;
      
      %% Reshape if its an array of objects
      if numel(decompOut)>1
        decompOut = reshape(decompOut,size(obj));
      end;
    end
    
    function isEmpty = isempty(obj)
      if isempty(obj.tfData)
        isEmpty = true;
      else
        isEmpty = false;
      end;
    end
    
    function out = plus(obj,b)      
      out = applyFcnToTFX(obj,b,@plus);      
      if isa(b,'tfDecomp')
      for i = 1:numel(out.chanLabels)
        if numel(b.chanLabels)==1, idxB = 1; else, idxB = i; end;
        out.chanLabels{i} = [obj.chanLabels{i} '+' b.chanLabels{idxB}];
      end
      end;
    end;
    
    function out = minus(obj,b)      
      out = applyFcnToTFX(obj,b,@minus);      
      if isa(b,'tfDecomp')
        for i = 1:numel(out.chanLabels)
          if numel(b.chanLabels)==1, idxB = 1; else, idxB = i; end;
          out.chanLabels{i} = [obj.chanLabels{i} '-' b.chanLabels{idxB}];
        end
      end;
    end;    
    
    function out = rdivide(obj,b)
      out = applyFcnToTFX(obj,b,@rdivide);   
      if isa(b,'tfDecomp')
      for i = 1:numel(out.chanLabels)
        if numel(b.chanLabels)==1, idxB = 1; else, idxB = i; end;
        out.chanLabels{i} = ['(' obj.chanLabels{i} ')/' b.chanLabels{idxB}];        
      end
      end;      
    end
    
    function out = times(obj,b)
      out = applyFcnToTFX(obj,b,@times);      
      if isa(b,'tfDecomp')
        for i = 1:numel(out.chanLabels)
          if numel(b.chanLabels)==1, idxB = 1; else, idxB = i; end;
          out.chanLabels{i} = [obj.chanLabels{i} '*' b.chanLabels{idxB}];
        end
      end;
    end;
  
    function out = power(obj,a)
      out = applyFcnToTFX(obj,a,@power);      
    end

    function out = cat(dim,obj,a,varargin)
      
      assert(isa(a,'tfDecomp'),'Can only concatenate with like objects');
      
      if isempty(obj)
        out = a;
        return;
      end;
      
      switch dim
        case 1
          error('Concatenation along frequency axis not implemented');
        case 2
          assert(isequal(obj.fVals,a.fVals),'Frequencies must match');
          assert(isequal(obj.chanLabels,a.chanLabels),'Labels must match');
          
          newType = ['TimeConcatenated'];
          
          concat = cat(2,obj.tfData,a.tfData);
          tVals = cat(1,obj.tVals,a.tVals); %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          
          chanLabels = obj.chanLabels; %#ok<PROPLC>
        case 3
          assert(isequal(obj.fVals,a.fVals),'Frequencies must match');
          assert(isequal(obj.tVals,a.tVals),'Times must match');

          newType = ['cat(' obj.decompType ',' a.decompType ')'];
          concat = cat(3,obj.tfData,a.tfData);          
          tVals = obj.tVals; %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          chanLabels = [obj.chanLabels ; a.chanLabels]; %#ok<PROPLC>
        otherwise
          error('Invalid dimension for concatenation');
      end;
       
      out = tfDecomp(newType,concat,tVals,fVals,chanLabels); %#ok<PROPLC>
      
      if ~isempty(varargin)
        out = cat(dim,out,varargin{:});
      end;
      
    end
    
    function out = sum(obj,dim)
      if ~exist('dim','var'), dim = 1; end;
      
      switch dim
        case 1          
          tVals = obj.tVals; %#ok<PROPLC>
          fVals = 0; %#ok<PROPLC>
          chanLabels = obj.chanLabels;           %#ok<PROPLC>
          newType = ['MeanF(' obj.decompType ')'];          
        case 2
          tVals = mean(obj.tVals); %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          chanLabels = obj.chanLabels; %#ok<PROPLC>
          newType = ['MeanT(' obj.decompType ')'];
        case 3
          tVals = obj.tVals; %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          chanLabels = {'MeanChan'}; %#ok<PROPLC>
          newType = ['MeanC(' obj.decompType ')'];
        otherwise
          error('Invalid dimension selection');
      end
      
      meanTFX = sum(obj.tfData,dim);
      out = tfDecomp(newType,meanTFX,tVals,fVals,chanLabels); %#ok<PROPLC>
      
    end
   function out = mean(obj,dim)
      if ~exist('dim','var'), dim = 1; end;
      
      switch dim
        case 1          
          tVals = obj.tVals; %#ok<PROPLC>
          fVals = 0; %#ok<PROPLC>
          chanLabels = obj.chanLabels;           %#ok<PROPLC>
          newType = ['MeanF(' obj.decompType ')'];          
        case 2
          tVals = mean(obj.tVals); %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          chanLabels = obj.chanLabels; %#ok<PROPLC>
          newType = ['MeanT(' obj.decompType ')'];
        case 3
          tVals = obj.tVals; %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          chanLabels = {'MeanChan'}; %#ok<PROPLC>
          newType = ['MeanC(' obj.decompType ')'];
        otherwise
          error('Invalid dimension selection');
      end
      
      meanTFX = mean(obj.tfData,dim);
      out = tfDecomp(newType,meanTFX,tVals,fVals,chanLabels); %#ok<PROPLC>
      
    end
    
    function out = std(obj,W,dim)
      if ~exist('dim','var'), dim = 1; end;
      if ~exist('W','var'), W = 0; end;
                  
      switch dim
        case 1          
          tVals = obj.tVals; %#ok<PROPLC>
          fVals = 0; %#ok<PROPLC>
          chanLabels = obj.chanLabels;           %#ok<PROPLC>
          newType = ['StdF(' obj.decompType ')'];          
        case 2
          tVals = 0; %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          chanLabels = obj.chanLabels; %#ok<PROPLC>
          newType = ['StdT(' obj.decompType ')'];
        case 3
          tVals = obj.tVals; %#ok<PROPLC>
          fVals = obj.fVals; %#ok<PROPLC>
          chanLabels = {'StdChan'}; %#ok<PROPLC>
          newType = ['StdC(' obj.decompType ')'];
        otherwise
          error('Invalid dimension selection');
      end      
      
      stdTFX = std(obj.tfData,W,dim);
      out = tfDecomp(newType,stdTFX,tVals,fVals,chanLabels); %#ok<PROPLC>
      
    end
    
    
  end
  
end

