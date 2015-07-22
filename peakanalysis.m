function [minTime, minValue, maxTime, maxValue] = peakanalysis(time, voltage, threshold)

% Find and analyze time and value of peaks in a trace
%
% Input:
%   time = Time of a peak
%   voltage = Voltage of the peak
%   threshold = Threshold above which a peak must be to be analyzed;
%               set this based on what it looks like it should be from
%               a preliminary plot of your data
%
%
% Output:
%   minTime = Vector of times at which minimum peaks occur
%   minValue = Vector of value at minimum peaks
%   maxTime = Vector of times at which maximum peaks occur
%   maxValue = Vector of values at maximum peaks

minTime = [];
minValue = [];
maxTime = [];
maxValue = [];

% Use peakdet function to find mins and maxs
[maxtab, mintab] = peakdet(voltage, threshold, time);

% Separate the mintab and maxtab into two respective
% arrays by variable (time of peak and value of peak)
if mintab ~= 0
    minTime = mintab(:,1);
    minValue = mintab(:,2);
elseif mintab == 0
    minTime = 0;
    minValue = 0;
end

if maxtab ~= 0
    maxTime = maxtab(:,1);
    maxValue = maxtab(:,2);
elseif maxtab == 0
    maxTime = 0;
    maxValue = 0;
end
    
% Plot the traces and distinctly label the maxima and minima
%plot(time,voltage)
%hold on
%plot(minTime,minValue,'bo')
%plot(maxTime,maxValue, 'ro')
%ylabel('Current (pA)')
%xlabel('Time')
%hold off


