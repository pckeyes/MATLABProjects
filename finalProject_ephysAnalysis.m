%% Chen Lab (First Rotation) Ephys Data Analysis
%  Created by: Piper Keyes
%  Created on: 11/21/14

%% Paired Pulse Facillitation Analysis Part 1
%Get the ppf ratio of each individual trace with:
% [rawData, samplingInterval, fileInfo] = abfload('FILENAME.abf');
smoothingFactor = 10; %Choose an appropriate smoothing factor for your traces
rawDataSmooth = moving_average(rawData(:,1,2),smoothingFactor); %Smooths your trace
figure(1)
plot(rawDataSmooth,'k') %Use this plot to files that are not too noisy and 
                        %don't have confounding spontaneous activity and 
                        %to choose a range of the data that only includes 
                        %stimulated responses
                         
%% Paired Pulse Facillitation Analysis Part 2
traceStart = 300;
traceLength = 5500;
dataRange = traceStart:traceLength; %Pick a range of data that only includes stim responses
threshold = 30; %Use an appropriate threshold to only choose the stim response peaks
[minTime,minValue,maxTime,maxValue] = ...
        peakanalysis(dataRange,rawDataSmooth(dataRange,1),threshold); 
     

figure(1)
plot(rawDataSmooth(dataRange,1),'k')
hold on
plot((minTime-traceStart),minValue,'mo') %Check for accurate minima here, and
                                   %change the subtracted number based on
                                   %your chosen range
minAmplitude = zeros(1,2);
cmap = ['r','b']; % to make each baseline have a different color
for iPeak = 1:length(minValue)
    baseline = mean(rawDataSmooth((minTime(iPeak)-125):(minTime(iPeak)-100),1));
    plot(linspace(0,traceLength,1000),baseline,'r','Color',cmap(iPeak))
    minAmplitude(iPeak) = minValue(iPeak) - baseline;
end
hold off
ppfTrace = minAmplitude(2)/minAmplitude(1);
%Any traces where the cell started bursting cannot be used. For useable
%cells, save these values in an excel file sub-divided into columns based 
%on inter-stimulation interval of that trace and take the average of the
%ratios. Add those values to a new array called ppfTotal which is also 
%subdivided into columns by cell and rows by interval.

%% Paired Pulse Facillitation Analysis Part 3

ppfTotal = [1.4771, 2.8593, 1.8988, 1.5402, 1.3613 %These data taken from excel
            1.3163, 2.3494, 1.3737, 1.2603, 1.3993
            0.9572, 1.4897, 1.2405, 1.0018, 1.0717
            1.0105, 1.1351, 1.0726, 0.8099, 0.9373];
        
intervalNumber = [1,2,3,4];
intervalName = [50,100,200,400];
amplitudeMean = mean(ppfTotal,2);

%  Calculate standard error
stddevAmplitude = std(ppfTotal,0,2);
stderrorAmplitude = zeros(1,4);
for iGroup = 1:4
    stderrorAmplitude(iGroup) = ...
        stddevAmplitude(iGroup)/sqrt(length(ppfTotal(iGroup,:)));
end

figure(2)
errorb(intervalNumber, amplitudeMean, stderrorAmplitude, 'Color', 'k')
hold on
plot(intervalNumber, amplitudeMean,'c.','MarkerSize',30)
set(gca,'XTickLabel', {intervalName},...
        'XTick',[1 2 3 4]);
xlabel('Inter-stimulation Interval (ms)')
ylabel('Peak 2/Peak 1 Ratio')
title('Paired Pulse Facilitation in CA1 Pyramidal Neurons')
hold off

%% Spontaneous Activity Part 1

% Load whatever file you want using this code:
% [rawData, samplingInterval, fileInfo] = abfload('FILENAME.abf');

% Separate big file into a matrix with all the traces and plot them
traceLength = 50000; %Length of the total trace or wherever you want analysis to end
traceStart = 1500; %Choose a starting point that will exclude the stimulation artifact
smoothingFactor = 10;
numTraces = 30; %Number of traces to be analyzed form the file
allTraces = zeros(traceLength,numTraces);
traceTime = 0:10; %For converting the data to seconds on graphs
close all
for iTrace = 1:numTraces
    rawDataSmooth = moving_average(rawData(:,1,iTrace),smoothingFactor); %Smooths your trace
    allTraces(:,iTrace) = rawDataSmooth;
    figure(iTrace)
    plot(allTraces(:,iTrace),'k','LineWidth',.05)
    xlim([traceStart,traceLength]) %This range excludes timepoints containing the 
                       %stimulation artifact
    ylim([-500,10])  %After running through this the first time, check this
                        %axis and choose a reasonable range
    xlabel('Time (s)')
    ylabel('Current (pA)')
    set(gca,'XTickLabel', (traceTime),...
        'XTick',linspace(0,traceLength,11));
    hold on
end

%% Spontaneous Activity Part 2

%Count and label the peaks (minima) in each trace and save the values in an
%excel file based on condition
dataRange = traceStart:traceLength; %Same range of data used in the xlim command above
threshold = 10; %Use an appropriate threshold to only choose true current peaks
thisMinCount = zeros(1,numTraces);
cellMinCount = 0;
for iTrace = 1:numTraces
    [minTime,minValue,maxTime,maxValue] = ...
        peakanalysis(dataRange,allTraces(dataRange,iTrace),threshold);
    if min(allTraces(dataRange,iTrace)) > -1000
        if length(minValue) > 0
            thisMinCount(iTrace) = length(minValue);
            cellMinCount = cellMinCount + thisMinCount(iTrace);
            figure(iTrace)
            plot(minTime,minValue,'mo')
            hold off
        elseif length(minValue) == 0
            cellMinCount = cellMinCount + 0;
        end
    elseif min(allTraces(dataRange,iTrace)) <= -1000 %Removes files with bursting from analysis
        cellMinCount = cellMinCount + 0;
    end
end

%Calculate average amplitude for each cell and save the values in an excel
%file based on condition.
thisAveAmplitude = zeros(1,numTraces);
tracesWithPeaks = 0;
cmap = hsv(20); % to make each baseline have a different color
for iTrace = 1:numTraces
            [minTime,minValue,maxTime,maxValue] = ...
        peakanalysis(dataRange,allTraces(dataRange,iTrace),threshold);
     if min(allTraces(dataRange,iTrace)) > -1000    
        if length(minValue) > 0
            for iPeak = 1:length(minValue)
            baseline = mean(allTraces((minTime(iPeak)-500):(minTime(iPeak)-50),iTrace),1);
            figure(iTrace)
            hold on
            %plot(linspace(traceStart,traceLength,5000),baseline,'r','Color',cmap(iPeak,:))   %Check baseline validity,
            %linspace dimensions based
            %on your data range
            hold off
            thisAveAmplitude(iTrace) = mean(minValue-baseline,1);
            end
        elseif length(minValue) == 0
            thisAveAmplitude(iTrace) = 0;
        end
     elseif min(allTraces(dataRange,iTrace)) <= -1000 %Removes files with bursting from analysis
        thisAveAmplitude(iTrace) = 0;
    end
end
for iTrace = 1:30
    if thisMinCount(iTrace) > 0
       tracesWithPeaks = tracesWithPeaks + 1;
    elseif thisMinCount(iTrace) == 0
       tracesWithPeaks = tracesWithPeaks + 0;
    end
end
cellAveAmplitude = (sum(thisAveAmplitude))/tracesWithPeaks;

%Save cellMinCount and cellAveAmplitude into an excel file  
%% Spontaneous Activity Part 3

%Plot peak counts for control vs experimental
peakCount = [333 246 96 60 130     %Data taken from excel file
             83 257 1081 582 852];

group = [1,2];
groupName = cell(1,2);
groupName{1} = 'control';
groupName{2} = 'picrotoxin';
peakMean = mean(peakCount,2);

%  Calculate standard error
stddevPeak = std(peakCount,0,2);
stderrorPeak = zeros(1,2);
for iGroup = 1:2
    stderrorPeak(iGroup) = ...
        stddevPeak(iGroup)/sqrt(length(peakCount(iGroup,:)));
end

%Create graph of means and run stats comparing the two
subplot(1,2,2)
bar(group, peakMean)
hold on
errorb(group, peakMean, stderrorPeak, 'Color', 'k')
set(gca,'XTickLabel', (groupName),...
        'XTick',[1 2]);
xlim([0,3])
xlabel('Drug Condition')
ylabel('# Spontaneous Currents')
title('Spontaneous Activity in CA1 Pyramidal Neurons')
hold off

[peakRejectNull,peakPValue,peakCI] = ttest2(peakCount(1,:),peakCount(2,:));


%Plot average amplitude for control vs experimental
aveAmplitude = [-10.2737, -10.8044, -10.813, -9.5073, -6.7164  %This data taken from excel file
                -5.4764, -6.7108, -16.8754, -15.9432, -29.3556];

group = [1,2];
groupName = cell(1,2);
groupName{1} = 'control';
groupName{2} = 'picrotoxin';
amplitudeMean = mean(aveAmplitude,2);

%  Calculate standard error
stddevAmplitude = std(aveAmplitude,0,2);
stderrorAmplitude = zeros(1,2);
for iGroup = 1:2
    stderrorAmplitude(iGroup) = ...
        stddevAmplitude(iGroup)/sqrt(length(aveAmplitude(iGroup,:)));
end

%Create graph of means and run stats comparing the two
subplot(1,2,1)
bar(group, amplitudeMean)
hold on
errorb(group, amplitudeMean, stderrorAmplitude, 'Color', 'k')
set(gca,'XTickLabel', (groupName),...
        'XTick',[1 2]);
xlim([0,3])
xlabel('Drug Condition')
ylabel('Amplitude (pA)')
title('Spontaneous Activity in CA1 Pyramidal Neurons')
hold off

[amplitudeRejectNull,amplitudePValue,amplitudeCI] = ttest2(aveAmplitude(1,:),aveAmplitude(2,:));