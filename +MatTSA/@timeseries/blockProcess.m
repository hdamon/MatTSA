function output = blockProcess(tSeriesIn,procHandle,varargin)
% Process a MatTSA.timeseries object in temporal blocks
%
% output = blockProcess(tSeriesIn,procHandle,varargin)
%
% Inputs
% ------
%   tSeriesIn : EEG to process
%  procHandle : Function handle to the analysis function
%                 Must take a single MatTSA.timeseries object as input, 
%                 and return a structure as the output. The fields of the
%                 structure should be the same for all input EEG's, so that
%                 the output blocks can be appropriately concatenated.
%
% Param-Value Inputs
% ------------------
%   'blockSize'    : Size of block to process (in samples)
%   'blockOverlap' : Amount of overlap in blocks (in samples)
%                     THIS OPTION NOT CURRENTLY IMPLEMENTED
%   'catDim'       : Dimension to concatenate outputs along;%
%                     DEFAULT: 2
%   'outputType'   : Type of output
%                     DEFAULT: 'MatTSA.timeseries'
%
% Output
% -------
%  output : Concatenation of processed blocks.
%
%

p = inputParser;
p.addParameter('blockSize',[],@(x) isnumeric(x)&&isscalar(x));
p.addParameter('blockOverlap',[],@(x) isnumeric(x)&&isscalar(x));
% Should probably remove these two options.
p.addParameter('catDim',2);
p.addParameter('outputType','MatTSA.tfDecomp', @(x) isempty(x)||ischar(x));
p.parse(varargin{:});

blockSize    = p.Results.blockSize*tSeriesIn.sampleRate;
blockOverlap = p.Results.blockOverlap;

nBlocks = floor(size(tSeriesIn,1)/blockSize);

output = eval(p.Results.outputType);

%% Loop Over Blocks
tic;
for idxBlock = 1:nBlocks
  currTime = toc;
  avgTime = currTime/idxBlock;
  compTime = avgTime*(nBlocks-idxBlock);
  dSize = fprintf(['Processing block ' num2str(idxBlock) ' of ' num2str(nBlocks) ' Total Time: ' num2str(currTime) ' Avg/Block: ' num2str(avgTime) ' Est Time Remaining: ' num2str(compTime)]);  
  offset = (idxBlock-1)*blockSize;
  
  % Get the Block
  s.type = '()';
  s.subs = {offset+1:offset+blockSize,':'};
  tmpEEG = subsref(tSeriesIn,s);
  
  % This Doesn't work because inside object functions, Matlab tries the
  % builtin function by default.
  % tmpEEG = tSeriesIn(offset+1:offset+blockSize,:);
  
  % Process the Block
  tmp = procHandle(tmpEEG);
  if ~isempty(tmp)
    % Discard failed blocks. Although procHandle should really be returning
    % a valid, but empty, block.
    output = cat(2,output,tmp);
  end;
  fprintf(repmat('\b',1,dSize));
end
disp('Completed block processing');

%% If the output of each block is a struct, concatenate the individual fields
if isstruct(output)
  f = fields(output);
  for idxF = 1:numel(f)
    tmp.(f{idxF}) = cat(2,output.(f{idxF}));
  end
  output = tmp;
end




































