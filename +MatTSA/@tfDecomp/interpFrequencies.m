function tfOut = interpFrequencies(tfIn,freqOut,type)

% Interpolate Frequencies in a MatTSA.tfDecomp object
%
% tfOut = interpFrequencies(tfIn,freqOut,type)
%
% Inputs
% ------
%     tfIn : Input MatTSA.tfDecomp object
%  freqOut : Output frequencies
%             Can be defined either as a vector of individual frequencies:
%             (ie: freqOut = [1:0.5:150])
%             or as a cell array of frequency bands to interpolate and
%             average/sum over:
%             (ie: freqOut = { [0.25:0.25:5], [5.25:0.25:10], .... }
%             The final output decomposition will have values at a total of
%             numel(freqOut) frequencies.
%     type : Post-processing type for cell-array inputs.
%             Options: {'sum', 'mean'}  Default: 'mean'
%             When the output frequencies are provided as a cell array, for
%             each cell, the input is interpolated to each value in the
%             cell, and then either averaged ('mean') or summed ('sum') to
%             produce the final output value.
%
% Outputs
% -------
%   tfOut : Output MatTSA.tfDecomp object
%
% WARNING: Interpolating and/or averaging complex valued tfDecomp objects
%             may give inconsistent results.
%
%

if ~exist('type','var'), type = 'mean'; end;

if isequal(tfIn.dataType,'complex')
%  warning('Interpolating and averaging complex time frequency values is not advised');
end

% The first two dimensions are frequency and time. Recurse over all other
% dimensions
indexOver = size(tfIn);
indexOver = indexOver(3:end);

% Strip trailing singleton dimensions
while ~isempty(indexOver)&&(indexOver(end)==1)
  indexOver = indexOver(1:end-1);
end

nIndex = numel(indexOver);

% Recurse over non-singleton dimensions
if nIndex>0
  
  allIndex(1:(nIndex+1)) = {':'};
  
  tfOut = MatTSA.tfDecomp;
  for idxLoop = 1:indexOver(nIndex)
    newIndex = [ allIndex {idxLoop}];
    s.type = '()';
    s.subs = newIndex;
    tfTmp = subsref(tfIn,s);
    tfTmp = interpFrequencies(tfTmp,freqOut,type);
    tfOut = cat(nIndex+2,tfOut,tfTmp);
  end
  return;
end % END non-singleton recursion

% At this point, tfIn.tfData should be a 2D matrix.

if iscell(freqOut)
  % Averaging over bands
  %disp('Interpolating and averaging in bands');
  Fin = tfIn.fVals;
  
  clear Fout
  data = tfIn.tfData; % to prevent multiple referencing calls
  Fout = nan(numel(freqOut),1);
  interpWAVE = nan(numel(freqOut),size(tfIn,2));
  for idxF = 1:numel(freqOut)
    %disp(['    Interpolating band #' num2str(idxF)]);
    Fout(idxF) = mean(freqOut{idxF});       
    switch type
      case 'mean'        
        interpWAVE(idxF,:) = mean(interp1(Fin,data,freqOut{idxF}),1);
      case 'sum'        
        interpWAVE(idxF,:) = sum(interp1(Fin,data,freqOut{idxF}),1);
    end
  end  
else
  % Interpolating to individual frequencies
  %disp('Interpolating to individual frequencies');
  Fin = tfIn.fVals; % Input Frequencies
  Fout = freqOut;
  
  interpWAVE = interp1(Fin,tfIn.tfData,Fout);    
end

interpParams.inputFreq = Fin;
interpParams.outputFreq = Fout;
interpParams.freqDef = freqOut;
decompParams = {tfIn.decompParams {interpParams}};

tfOut = MatTSA.tfDecomp(interpWAVE,'decompType',tfIn.decompType,...
  'dataType',tfIn.dataType,...
  'tVals',tfIn.tVals,...
  'fVals', Fout,...
  'chanLabels',tfIn.chanLabels,...
  'decompParams',decompParams);

end % END interpFrequencies




