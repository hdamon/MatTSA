function obj = subsasgn(obj,s,varargin)
% Overloaded subsasgn for crlBase.type.timeseries objects
%
   switch s(1).type
      case '.'
         if length(s) == 1
            % Implement obj.PropertyName = varargin{:};
            obj = builtin('subsasgn',obj,s,varargin{:});
         elseif length(s) == 2 && strcmp(s(2).type,'()')
            % Implement obj.PropertyName(indices) = varargin{:};
            if (numel(obj)==1) && (isequal(s(1).subs,'data'))
              if numel(s(2).subs)==2
                rowIdx = s(2).subs{1};
                colIdx = getElectrodeIndex(s(2).subs{2},true);
                obj.array_(rowIdx,colIdx) = varargin{:};
              else
              keyboard;
              end;
            else
              obj = builtin('subsasgn',obj,s,varargin{:});
            end;
         else
            % Call built-in for any other case
            obj = builtin('subsasgn',obj,s,varargin{:});
         end
      case '()'
         if length(s) == 1
            % Implement obj(indices) = varargin{:};
            obj = builtin('subsasgn',obj,s,varargin{:});
         elseif length(s) == 2 && strcmp(s(2).type,'.')
            % Implement obj(indices).PropertyName = varargin{:}{:};
            obj = builtin('subsasgn',obj,s,varargin{:});
         elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
            % Implement obj(indices).PropertyName(indices) = varargin{:};
            if numel(obj)>1
              tmp = subsref(obj,s(1));
              obj(s(1).subs{:}) = subsasgn(tmp,s(2:end),varargin{:});
            else
              obj = builtin('subsasgn',obj,s,varargin{:});
            end;
         else
            % Use built-in for any other expression
            obj = builtin('subsasgn',obj,s,varargin{:});
         end       
      case '{}'
         if length(s) == 1
            % Implement obj{indices} = varargin{:}
            obj = builtin('subsasgn',obj,s,varargin{:});
         elseif length(s) == 2 && strcmp(s(2).type,'.')
            % Implement obj{indices}.PropertyName = varargin{:}
            obj = builtin('subsasgn',obj,s,varargin{:});
            % Use built-in for any other expression
            obj = builtin('subsasgn',obj,s,varargin{:});
         end
      otherwise
         error('Not a valid indexing expression')
   end


  function outIdx = getElectrodeIndex(valIn,isNumericValid)
    if isequal(valIn,':')
      outIdx = ':';
    elseif iscellstr(valIn) % Cell array of strings
      outIdx = getIdxFromStringCell(valIn);
    elseif ischar(valIn) % Single character string
      outIdx = getIdxFromStringCell({valIn});
    else
      % Poorly defined: Can't concreately determine if its supposed
      % to index into time or channels
      if isNumericValid
        outIdx = valIn;
      else
        error('Invalid Numeric Indexing');
      end;
    end;
  end

  function outIdx = getIdxFromStringCell(cellIn)
    %%
    outIdx = zeros(1,numel(cellIn));
    for idx = 1:numel(outIdx)
      tmp = find(strcmp(cellIn{idx},obj.chanLabels));
      assert(numel(tmp)==1,'Invalid Index Labels');
      outIdx(idx) = tmp;
    end
  end

end