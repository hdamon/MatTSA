function n = numArgumentsFromSubscript(obj,s,indexingContext)
  

n = numArgumentsFromSubscript@labelledArray(obj,s,indexingContext);
return;
% Default
if ~iscell(s(1).subs)
  n = builtin('numArgumentsFromSubscript',obj,s,indexingContext);
else
  n = 1;
end;
  
% % Optional overrides
%    switch indexingContext
%       case matlab.mixin.util.IndexingContext.Statement
%         switch s(1).type
%           case '()'
%             if numel(obj)==1
%              % When 
%              n = 1; % nargout for indexed reference used as statement            
%             end  
%           otherwise
%         end;            
%       case matlab.mixin.util.IndexingContext.Expression
%         switch s(1).type
%           case '()'
%             if numel(obj)==1
%              % When 
%              n = 1; % nargout for indexed reference used as statement            
%             end  
%           otherwise
%         end;         
%          %n = 1; % nargout for indexed reference used as function argument
%       case matlab.mixin.util.IndexingContext.Assignment        
%         
%          %n = 1; % nargin for indexed assignment
%    end
% end