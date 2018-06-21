function varargout = subsref(obj,s)
% subsref method for MatTSA.timeseries
%
%
% There have been significant modifications to the way the timeseries
% object uses indexing.
%

switch s(1).type
  case '.'
    %%
    if (length(s)==2)&isequal(s(2).type,'()')     
            
      if (numel(obj)==1)&(isequal(s(1).subs,'tfX'))        
        % Enables use of non-numeric referencing for obj.data(a,b) type
        % referencing.
        %
        % IE:
        %  obj.data('Cz')
        %  obj.data(1:10,{'Cz' 'Pz'});
        tmp = obj.subsref(s(2));
        varargout = {tmp.tfX};
      else
        if nargout==0
          builtin('subsref',obj,s);
        else
          varargout = {builtin('subsref',obj,s)};
        end;
      end;
    else
      varargout = {builtin('subsref',obj,s)};
    end;
  case '()'
    %% Parenthetical Referencing
    if length(s) == 1
      %% Implement obj(indices)
      if numel(obj)==1
        % Internal object indices can only be accessed individually
        if numel(s.subs)==1
          %% To avoid ambiguity, single indexing into timeseries must be 
          %  with a character string. 
          
          if isequal(s.subs{1},1)
            % Permits referencing foo(1) when foo is a single
            % tfDecomp object
            varargout = {obj};
            return;
          end;
          
          fIdx = ':';
          tIdx = ':';
          chanIdx = crlEEG.util.getIndexIntoCellStr(obj.labels,s.subs{1},false);
        elseif numel(s.subs)==2
          if numel(size(obj.tfX))>2
            error('Poorly defined indexing expression');          
          else
            fIdx = s.subs{1};
            tIdx = s.subs{2};
            chanIdx = 1;
          end;
        elseif numel(s.subs)==3
          fIdx = s.subs{1};
          tIdx = s.subs{2};
          chanIdx = crlEEG.util.getDimensionIndex(obj.labels,s.subs{3},true);          
        end;
        
        if islogical(fIdx) && (numel(fIdx)==size(obj,1))
          fIdx = find(fIdx);
        end;
        
        if islogical(tIdx) && (numel(tIdx)==size(obj,2))
          tIdx = find(tIdx);
        end;
        
        % Time and Frequency References Must be Numeric At this Point
        assert(isequal(fIdx,':')||isnumeric(fIdx),...
          'First dimension indexing must be numeric');
        
        assert(isequal(tIdx,':')||isnumeric(tIdx),...
          'Second dimension indexing must be numeric');
        
        % Output a new timeseries object
        varargout = {obj.subcopy(fIdx,tIdx,chanIdx)};
      else
        % Indexing into an array of tfDecomp objects
        % returns a tfDecomp object.
        varargout = {obj(s.subs{:})};
      end;
    elseif length(s) == 2 && strcmp(s(2).type,'.')
      %% Implement obj(ind).PropertyName
      if ( numel(s(1).subs)==1 || numel(s(1).subs)==2 )
        tmp = obj.subsref(s(1));
        varargout = {builtin('subsref',tmp,s(2))};
      else
        varargout = {builtin('subsref',obj,s)};
      end;
    elseif ( length(s) == 3 && strcmp(s(2).type,'.') ...
                                && strcmp(s(3).type,'()') )
      %% Implement obj(indices).PropertyName(indices)
      tmp = subsref(obj,s(1));      
      varargout = {subsref(tmp,s(2:3))};
    else
      %% Use built-in for any other expression
      varargout = {builtin('subsref',obj,s)};
    end
  case '{}'
    %%
    varargout = {builtin('subsref',obj,s)};
    
  otherwise
    error('Not a valid indexing expression')
end


end


