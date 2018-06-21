function outTseries = filter(tseries,dFilter)
% Overloaded filter function for MatTSA.timeseries objects
%
% Inputs
% ------
%  tseries : A MatTSA.timeseries object to be filtered
% dFilter : A matlab digitalFilter object (typically created with
% designfilt)
%

tmp = filter(dFilter,tseries.data(:,tseries.getChannelsByType('data')));

outTseries = tseries.copy;
outTseries.data(:,outTseries.getChannelsByType('data')) = tmp;

end