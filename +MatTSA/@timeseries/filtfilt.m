function outTseries = filtfilt(tseries,dFilter)
% Overloaded filtfilt function for crltseries.type.timeseries objects
%
% Inputs
% ------
%   tseries : A crltseries.type.tseries.object to be filtered
%  dFilter : A Matlab digital filter (typically created with designfilt)
%



dataChans = tseries.getChannelsByType('data');

tmp = filtfilt(dFilter,tseries.data(:,dataChans));

outTseries = tseries.copy;
outTseries.data(:,dataChans) = tmp;

end
