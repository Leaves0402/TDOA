%% TDOA Localization Error Heatmap - Dynamic Area Simulation
% 本版本已整合使用者自訂定位面積、Linear Least Squares 初始定位、
% 虛擬反矩陣（pinv）遠場降秩修正，以及系統取樣率（PRR）量化誤差模擬。

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

% 1. 使用者輸入參數
area_side_m = input('請輸入定位範圍總邊長 (單位:公尺, 例如 100 或 200): ');
grid_spacing_m = input('請輸入格線間距 (單位:公尺, 例如 5): ');
PRR_Hz = input('請輸入系統取樣率 PRR (單位:Hz, 例如 20000): ');

if isempty(area_side_m) || area_side_m <= 0
    area_side_m = 200; % 預設值
end
half_side = area_side_m / 2;

if grid_spacing_m <= 0 || PRR_Hz <= 0
    error('數值必須大於 0');
end

% 根據輸入的面積動態建立邊界與中心點網格
xEdges = -half_side:grid_spacing_m:half_side;
yEdges = -half_side:grid_spacing_m:half_side;

xCenters = xEdges(1:end-1) + grid_spacing_m/2;
yCenters = yEdges(1:end-1) + grid_spacing_m/2;

%% Build error maps
[errorPercentMap, estimatedXMap, estimatedYMap] = buildErrorMap( ...
    sensors, c, xCenters, yCenters, refIndex, PRR_Hz);

classMap = classifyErrorMap(errorPercentMap);

%% Plot and save heatmap
fig = plotErrorHeatmap(classMap, xCenters, yCenters, sensors, grid_spacing_m, half_side);
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

%% 線性 Least Squares 計算函數 (第一階段)
function p0 = computeInitialPositionLLS(sensors, measuredTDOA, c, refIndex)
    numSensors = size(sensors, 1);
    s_ref = sensors(refIndex, :);
    
    A = zeros(numSensors - 1, 3);
    b = zeros(numSensors - 1, 1);
    
    idx = 1;
    for i = 1:numSensors
        if i == refIndex
            continue;
        end
        s_i = sensors(i, :);
        rho_i = c * measuredTDOA(idx); 
        
        A(idx, 1) = 2 * (s_ref(1) - s_i(1));
        A(idx, 2) = 2 * (s_ref(2) - s_i(2));
        A(idx, 3) = -2 * rho_i;
        
        b(idx, 1) = s_ref(1)^2 + s_ref(2)^2 - s_i(1)^2 - s_i(2)^2 + rho_i^2;
        
        idx = idx + 1;
    end
    
    % 使用虛擬反矩陣 (pinv) 穩定處理遠場造成的降秩/近奇異矩陣警告
    u = pinv(A) * b;
    
    % 取得線性代數初始猜測點 p0 = [x, y]
    p0 = [u(1), u(2)];
end 

function [estimatedPosition, exitFlag] = estimatePosition( ...
    measuredTDOA, sensors, c, refIndex, initialGuess)

    % 優先使用 lsqnonlin (Levenberg-Marquardt)
    if exist('lsqnonlin', 'file') == 2
        options = optimoptions('lsqnonlin', ...
            'Algorithm', 'levenberg-marquardt', ...
            'Display', 'off', ...
            'MaxFunctionEvaluations', 5000, ...
            'MaxIterations', 2000, ...
            'FunctionTolerance', 1e-12, ...
            'StepTolerance', 1e-12);

        [estimatedPosition, ~, ~, exitFlag] = lsqnonlin( ...
            @(p) tdoaResiduals(p, measuredTDOA, sensors, c, refIndex), ...
            initialGuess, [], [], options);
    else
        % 備案：如果沒有安裝最佳化工具箱，降級使用 fminsearch
        options = optimset( ...
            'Display', 'off', ...
            'MaxIter', 2000, ...
            'MaxFunEvals', 5000, ...
            'TolX', 1e-12, ...
            'TolFun', 1e-12);

        objectiveFunction = @(p) sum(tdoaResiduals(p, measuredTDOA, sensors, c, refIndex).^2);

        [estimatedPosition, ~, exitFlag] = fminsearch( ...
            objectiveFunction, initialGuess, options);
    end
end

function [errorPercentMap, estimatedXMap, estimatedYMap] = buildErrorMap( ...
    sensors, c, xCenters, yCenters, refIndex, PRR_Hz)

    numX = numel(xCenters);
    numY = numel(yCenters);

    errorPercentMap = NaN(numY, numX);
    estimatedXMap = NaN(numY, numX);
    estimatedYMap = NaN(numY, numX);

    for iy = 1:numY
        for ix = 1:numX
            sourceTrue = [xCenters(ix), yCenters(iy)];

            % 1. 計算真實到達時間 (TOA)
            trueArrivalTimes = computeArrivalTimes(sourceTrue, sensors, c);
            
            % 2. 模擬硬體取樣率造成的量化誤差 (Quantization)
            quantizedArrivalTimes = round(trueArrivalTimes * PRR_Hz) / PRR_Hz;
            
            % 3. 計算量化後的量測 TDOA
            sensorIndices = setdiff(1:size(sensors, 1), refIndex);
            measuredTDOA = quantizedArrivalTimes(sensorIndices) - quantizedArrivalTimes(refIndex);

            % 4. 進行位置估測 (演算法兩階段串連)
            % 階段一：使用 Linear Least Squares 求得初始位置 p0
            p0 = computeInitialPositionLLS(sensors, measuredTDOA, c, refIndex);
            
            % 若線性解發散或產生 NaN，則跳過此點
            if any(~isfinite(p0))
                continue;
            end

            % 階段二：將 p0 作為合理起始點，投入 LM 演算法進行非線性最佳化修正
            [estimatedPosition, exitFlag] = estimatePosition( ...
                measuredTDOA, sensors, c, refIndex, p0);
            
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

function fig = plotErrorHeatmap(classMap, xCenters, yCenters, sensors, grid_spacing_m, half_side)
    fig = figure('Color', 'w');

    imageHandle = imagesc(xCenters, yCenters, classMap);
    set(imageHandle, 'AlphaData', ~isnan(classMap));
    set(gca, 'YDir', 'normal');
    axis equal;
    axis tight;
    grid on;
    ax = gca;
    
    % 動態計算主要刻度步長（將畫面均勻分成 4 個區間，呈現 5 個主要標籤數字）
    tick_step = half_side / 2;
    ax.XTick = -half_side:tick_step:half_side;
    ax.YTick = -half_side:tick_step:half_side;

    % 開啟次要網格線，並動態設定為符合輸入的格線間距
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.XAxis.MinorTickValues = -half_side:grid_spacing_m:half_side;
    ax.YAxis.MinorTickValues = -half_side:grid_spacing_m:half_side;
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
    cb.Label.String = 'Localization Error (Percent of True Distance)';

    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title(sprintf('TDOA Error Heatmap (Area: %dm x %dm, Spacing: %.1fm)', ...
        half_side*2, half_side*2, grid_spacing_m), 'FontSize', 11, 'FontWeight', 'bold');

    hold on;
    % 繪製感測器位置
    plot(sensors(:, 1), sensors(:, 2), 'ko', ...
        'MarkerFaceColor', 'w', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    hold off;
end