function outTseries = filtfilt(tseries,dFilter)
% Overloaded filtfilt function for crltseries.type.timeseries objects
%
% Inputs
% ------
%   tseries : A crltseries.type.tseries.object to be filtered
%  dFilter : A Matlab digital filter (typically created with designfilt)
%


tmp = filtfilt(dFilter,tseries.data(:,tseries.getChannelsByType('data')));

outTseries = tseries.copy;
outTseries.data(:,outTseries.getChannelsByType('data')) = tmp;

end