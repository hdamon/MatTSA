function varargout = subsref(obj,s)
% subsref method for timeseries objects
%

switch s(1).type
  case '.'    
    
    if isequal(s(1).subs,'data')&&(numel(s)>1)&&isequal(s(2).type,'()')
      % Permits obj.data(<ref>) indexing.
      tmp = subsref(obj,s(1));
      s(2).subs = obj.getNumericIndex(s(2).subs{:});
      %varargout = {subsref(tmp,s(2:end))};      
      [varargout{1:nargout}] = subsref(tmp,s(2:end));
    elseif ismember(s(1).subs,{'addChannel' 'removeChannel'})
      % Bit of a hack. Is this really necessary? I think it should be
      % caught by the methods block of labelledArray.subsref.
      builtin('subsref',obj,s);
      %varargout = {obj};
      [varargout{1:nargout}] = obj;
    else
      %varargout = {subsref@labelledArray(obj,s)};
      [varargout{1:nargout}] = subsref@labelledArray(obj,s);
    end
           
  case '()'
    %varargout = {subsref@labelledArray(obj,s)};
    [varargout{1:nargout}] = subsref@labelledArray(obj,s);
    
  case '{}'
    %varargout = {subsref@labelledArray(obj,s)};
    [varargout{1:nargout}] = subsref@labelledArray(obj,s);
    
  otherwise
    error('Not a valid indexing expression')
end
