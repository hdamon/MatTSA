function varargout = subsref(obj,s)
% subsref method for timeseries objects
%

switch s(1).type
  case '.'    
    
    if isequal(s(1).subs,'tfData')&&(numel(s)>1)&&isequal(s(2).type,'()')
      % Permits obj.data(<ref>) indexing.
      tmp = subsref(obj,s(1));
      s(2).subs = obj.getNumericIndex(s(2).subs{:});
      varargout = {subsref(tmp,s(2:end))};      
    else
      varargout = {subsref@labelledArray(obj,s)};
    end
           
  case '()'
    varargout = {subsref@labelledArray(obj,s)};
    
  case '{}'
    varargout = {subsref@labelledArray(obj,s)};
    
  otherwise
    error('Not a valid indexing expression')
end
