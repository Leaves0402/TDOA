%% TDOA Localization Error Heatmap - Initial 200 m Simulation
% This initial version uses ideal TDOA values only. It does not include
% sampling-rate limits, time quantization error, or measurement noise.
%
% The initial guess is fixed at [1.0, 1.0] for every grid cell, so some
% source locations may converge to an incorrect solution.
%
% The sensor spacing is only 0.5 m while the localization area is
% 200 m x 200 m. Far-field regions may be unstable or have multiple
% feasible solutions.
%
% This is an initial version intended to be extended later into batch
% simulations for different sensor geometries.

clear; clc; close all;

%% Sensor and simulation settings
sensors = [
     0.0,  0.0;   % S1, reference sensor
     0.5,  0.0;   % S2
    -0.5,  0.0;   % S3
     0.0,  0.5;   % S4
     0.0, -0.5    % S5
];

refIndex = 1;
c = 343.0;        % Speed of sound, m/s

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptDir = pwd;
else
    scriptDir = fileparts(scriptPath);
end
repoRoot = fileparts(scriptDir);
resultsDir = fullfile(repoRoot, 'results');

if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

pngOutputPath = fullfile(resultsDir, 'tdoa_error_heatmap_initial_200m.png');
csvOutputPath = fullfile(resultsDir, 'tdoa_error_percent_map_initial_200m.csv');

xEdges = -100:5:100;
yEdges = -100:5:100;

xCenters = xEdges(1:end-1) + 2.5;
yCenters = yEdges(1:end-1) + 2.5;

initialGuess = [1.0, 1.0];

options = optimset( ...
    'Display', 'off', ...
    'MaxIter', 2000, ...
    'MaxFunEvals', 5000, ...
    'TolX', 1e-10, ...
    'TolFun', 1e-20);

%% Build error maps
[errorPercentMap, estimatedXMap, estimatedYMap] = buildErrorMap( ...
    sensors, c, xCenters, yCenters, refIndex, initialGuess, options);

classMap = classifyErrorMap(errorPercentMap);

%% Plot and save heatmap
fig = plotErrorHeatmap(classMap, xCenters, yCenters, sensors);
print(fig, pngOutputPath, '-dpng', '-r300');

%% Save numerical result
csvwrite(csvOutputPath, errorPercentMap);

%% Print statistics
validErrors = errorPercentMap(isfinite(errorPercentMap));

if isempty(validErrors)
    fprintf('minimum error percent: NaN\n');
    fprintf('maximum error percent: NaN\n');
    fprintf('mean error percent: NaN\n');
    fprintf('median error percent: NaN\n');
    fprintf('Maximum error location x, y: NaN, NaN\n');
    fprintf('Estimated position at maximum error: NaN, NaN\n');
    fprintf('True position at maximum error: NaN, NaN\n');
else
    minErrorPercent = min(validErrors);
    maxErrorPercent = max(validErrors);
    meanErrorPercent = mean(validErrors);
    medianErrorPercent = median(validErrors);

    finiteErrorMap = errorPercentMap;
    finiteErrorMap(~isfinite(finiteErrorMap)) = -Inf;
    [~, maxLinearIndex] = max(finiteErrorMap(:));
    [maxRow, maxCol] = ind2sub(size(errorPercentMap), maxLinearIndex);

    maxTruePosition = [xCenters(maxCol), yCenters(maxRow)];
    maxEstimatedPosition = [estimatedXMap(maxRow, maxCol), estimatedYMap(maxRow, maxCol)];

    fprintf('minimum error percent: %.12g\n', minErrorPercent);
    fprintf('maximum error percent: %.12g\n', maxErrorPercent);
    fprintf('mean error percent: %.12g\n', meanErrorPercent);
    fprintf('median error percent: %.12g\n', medianErrorPercent);
    fprintf('Maximum error location x, y: %.12g, %.12g\n', ...
        maxTruePosition(1), maxTruePosition(2));
    fprintf('Estimated position at maximum error: %.12g, %.12g\n', ...
        maxEstimatedPosition(1), maxEstimatedPosition(2));
    fprintf('True position at maximum error: %.12g, %.12g\n', ...
        maxTruePosition(1), maxTruePosition(2));
end

disp('Simulation completed.');

%% Local functions
function arrivalTimes = computeArrivalTimes(source, sensors, c)
    distances = sqrt(sum((sensors - source).^2, 2));
    arrivalTimes = distances ./ c;
end

function measuredTDOA = computeTDOA(source, sensors, c, refIndex)
    arrivalTimes = computeArrivalTimes(source, sensors, c);
    sensorIndices = setdiff(1:size(sensors, 1), refIndex);
    measuredTDOA = arrivalTimes(sensorIndices) - arrivalTimes(refIndex);
end

function residuals = tdoaResiduals(p, measuredTDOA, sensors, c, refIndex)
    modelTDOA = computeTDOA(p, sensors, c, refIndex);
    residuals = modelTDOA - measuredTDOA;
end

function [estimatedPosition, finalCost, exitFlag] = estimatePosition( ...
    measuredTDOA, sensors, c, refIndex, initialGuess, options)

    objectiveFunction = @(p) sum(tdoaResiduals( ...
        p, measuredTDOA, sensors, c, refIndex).^2);

    try
        [estimatedPosition, finalCost, exitFlag] = fminsearch( ...
            objectiveFunction, initialGuess, options);
    catch
        estimatedPosition = [NaN, NaN];
        finalCost = NaN;
        exitFlag = -1;
    end

    if numel(estimatedPosition) ~= 2 || any(~isfinite(estimatedPosition)) || ...
            ~isfinite(finalCost)
        estimatedPosition = [NaN, NaN];
        finalCost = NaN;
        exitFlag = -1;
    end
end

function [errorPercentMap, estimatedXMap, estimatedYMap] = buildErrorMap( ...
    sensors, c, xCenters, yCenters, refIndex, initialGuess, options)

    numX = numel(xCenters);
    numY = numel(yCenters);

    errorPercentMap = NaN(numY, numX);
    estimatedXMap = NaN(numY, numX);
    estimatedYMap = NaN(numY, numX);

    for iy = 1:numY
        for ix = 1:numX
            sourceTrue = [xCenters(ix), yCenters(iy)];

            measuredTDOA = computeTDOA(sourceTrue, sensors, c, refIndex);

            [estimatedPosition, ~, exitFlag] = estimatePosition( ...
                measuredTDOA, sensors, c, refIndex, initialGuess, options);

            if exitFlag <= 0 || any(~isfinite(estimatedPosition))
                continue;
            end

            xEst = estimatedPosition(1);
            yEst = estimatedPosition(2);
            xTrue = sourceTrue(1);
            yTrue = sourceTrue(2);

            errorM = sqrt((xEst - xTrue)^2 + (yEst - yTrue)^2);
            rTrue = sqrt(xTrue^2 + yTrue^2);
            errorPercent = errorM / max(rTrue, 1e-6) * 100;

            if isfinite(errorPercent)
                errorPercentMap(iy, ix) = errorPercent;
                estimatedXMap(iy, ix) = xEst;
                estimatedYMap(iy, ix) = yEst;
            end
        end
    end
end

function classMap = classifyErrorMap(errorPercentMap)
    classMap = NaN(size(errorPercentMap));

    classMap(errorPercentMap >= 0 & errorPercentMap < 10) = 1;
    classMap(errorPercentMap >= 10 & errorPercentMap < 100) = 2;
    classMap(errorPercentMap >= 100 & errorPercentMap < 1000) = 3;
    classMap(errorPercentMap >= 1000 & errorPercentMap < 10000) = 4;
    classMap(errorPercentMap >= 10000) = 5;
end

function fig = plotErrorHeatmap(classMap, xCenters, yCenters, sensors)
    fig = figure('Color', 'w');

    imageHandle = imagesc(xCenters, yCenters, classMap);
    set(imageHandle, 'AlphaData', ~isnan(classMap));
    set(gca, 'YDir', 'normal');
    axis equal;
    axis tight;
    grid on;
    box on;

    colormap([
        0.00, 0.45, 0.74;   % 0-10%
        0.47, 0.67, 0.19;   % 10-100%
        0.93, 0.69, 0.13;   % 100-1000%
        0.85, 0.33, 0.10;   % 1000-10000%
        0.49, 0.18, 0.56    % >10000%
    ]);

    caxis([0.5, 5.5]);

    cb = colorbar;
    cb.Ticks = 1:5;
    labelDash = char(8211);
    cb.TickLabels = { ...
        ['0' labelDash '10%'], ...
        ['10' labelDash '100%'], ...
        ['100' labelDash '1000%'], ...
        ['1000' labelDash '10000%'], ...
        '>10000%'};
    cb.Label.String = 'Localization error percent class';

    xlabel('x (m)');
    ylabel('y (m)');
    title('TDOA Localization Error Heatmap (5 Sensors, spacing = 0.5 m)');

    hold on;
    plot(sensors(:, 1), sensors(:, 2), 'ko', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 6);

    for iSensor = 1:size(sensors, 1)
        text(sensors(iSensor, 1) + 2.0, sensors(iSensor, 2) + 2.0, ...
            sprintf('S%d', iSensor), ...
            'Color', 'k', ...
            'FontWeight', 'bold', ...
            'BackgroundColor', 'w', ...
            'Margin', 1);
    end
    hold off;
end
