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
    if (length(s)==2)&&isequal(s(2).type,'()')
      
      if (numel(obj)==1)&&(isequal(s(1).subs,'data'))
        % Enables use of non-numeric referencing for obj.data(a,b) type
        % referencing.
        %
        % IE:
        %  obj.data('Cz')
        %  obj.data(1:10,{'Cz' 'Pz'});
        tmp = obj.subsref(s(2)); %This is lazy and slow when there are decompositions involved
        varargout = {tmp.data};
      else
        % This is poorly coded.
        if nargout==0
          % This in particular is really bad!
          builtin('subsref',obj,s);
        else
          %foo(1:nargout) = {builtin('subsref',obj,s)};
          varargout{1:nargout} = builtin('subsref',obj,s);
        end;
      end;
    elseif length(s)>1
      tmpObj = subsref(obj,s(1));
      varargout = {subsref(tmpObj,s(2:end))};
    else      
      varargout = {builtin('subsref',obj,s)};
    end;
  case '()'
    %%
    if length(s) == 1
      %% Implement obj(indices)
      if numel(obj)==1
        % Internal object indices can only be accessed individually
        if numel(s.subs)==1
          %% To avoid ambiguity, single indexing into timeseries must be
          %  with a character string.
          
          if isequal(s.subs{1},1)
            varargout = {obj};
            return;
          end;
          
          rowIdx = ':';
          colIdx = crlBase.util.getIndexIntoCellStr(obj.chanLabels,s.subs{1},false);
        elseif numel(s.subs)==2
          rowIdx = s.subs{1};
          colIdx = crlBase.util.getIndexIntoCellStr(obj.chanLabels,s.subs{2},true);
        else
          error('Invalid indexing expression');
        end;
        
        if (numel(rowIdx)==size(obj,1))
          rowIdx = logical(rowIdx);
          rowIdx = find(rowIdx);
        end;
        
        % Time must be indexed numerically
        assert(isequal(rowIdx,':')||isnumeric(rowIdx),...
          'Second dimension indexing must be numeric');
        
        % Output a new timeseries object
        varargout = {obj.subcopy(rowIdx,colIdx)};
      else
        % Indexing into individual timeseries must be done with individual
        % objects, and not arrays.
        varargout = {obj(s.subs{:})};
      end;
    elseif length(s) == 2 && strcmp(s(2).type,'.')
      % Implement obj(ind).PropertyName
      if ( numel(s(1).subs)==1 || numel(s(1).subs)==2 )
        tmp = obj.subsref(s(1));
        varargout = {builtin('subsref',tmp,s(2))};
      else
        varargout = {builtin('subsref',obj,s)};
      end;
    elseif ( length(s) == 3 && strcmp(s(2).type,'.') ...
        && strcmp(s(3).type,'()') )
      % Implement obj(indices).PropertyName(indices)
      if (numel(obj)>1) && strcmp(s(1).type,'()')
        tmp = obj.subsref(s(1));
        varargout = {tmp.subsref(s(2:end))};
      else        
        varargout = {builtin('subsref',obj,s)};
      end;
    else
      % Use built-in for any other expression
      varargout = {builtin('subsref',obj,s)};
    end
  case '{}'
    %%
    varargout = {builtin('subsref',obj,s)};
    
  otherwise
    error('Not a valid indexing expression')
end

end


