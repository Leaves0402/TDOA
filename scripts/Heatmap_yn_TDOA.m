%% 2D TDOA 定位誤差百分比熱力圖 (yn_TDOA 核心 + PRR 量化 + 新五大級距)
% 架構特色：
% 1. 演算法核心：採用 yn_TDOA 邏輯（無 LLS 線性預估、初始點固定為感測器幾何質心 p0）。
% 2. 空間維度：2D 平面，5 顆探頭呈十字形分佈（中央 S1 置於原點，其餘向外延伸 0.5m）。
% 3. 硬體限制：加入使用者自訂系統取樣率 (PRR_Hz)，模擬真實硬體的時間量化誤差。
% 4. 視覺化級距：全新五大誤差百分比分級 (<1%, 1-5%, 5-10%, 10-20%, >20%)。

clear; clc; close all;

%% ============================================================
% 1. 感測器配置與環境參數設定 (單位：公尺)
% ============================================================
sensors = [
     0.0,  0.0;   % S1, 參考感測器 (置於原點)
     0.5,  0.0;   % S2 (沿 +X 軸延展 0.5m)
    -0.5,  0.0;   % S3 (沿 -X 軸延展 0.5m)
     0.0,  0.5;   % S4 (沿 +Y 軸延展 0.5m)
     0.0, -0.5    % S5 (沿 -Y 軸延展 0.5m)
];

refIndex = 1;
c = 343.0;        % 聲速 (m/s)

% 透過 mean 計算出固定的初始猜測點 p0，稍微偏移 1 毫米 (0.001m)，避開原點奇異點
p0 = mean(sensors, 1)' + [0.001; 0.001];

%% ============================================================
% 2. 使用者互動輸入與網格建立
% ============================================================
area_side_m = input('請輸入定位範圍總邊長 (單位:公尺, 例如 100 或 200): ');
grid_spacing_m = input('請輸入格線間距 (單位:公尺, 例如 1 或 5): ');
PRR_Hz = input('請輸入系統取樣率 PRR (單位:Hz, 例如 20000 或 100000): ');

if isempty(area_side_m) || area_side_m <= 0, area_side_m = 200; end
if isempty(grid_spacing_m) || grid_spacing_m <= 0, grid_spacing_m = 5; end
if isempty(PRR_Hz) || PRR_Hz <= 0, PRR_Hz = 20000; end

half_side = area_side_m / 2;
xCenters = -half_side:grid_spacing_m:half_side;
yCenters = -half_side:grid_spacing_m:half_side;
[X_mesh, Y_mesh] = meshgrid(xCenters, yCenters);

errorPercentMap = zeros(size(X_mesh)); % 儲存每個網格點的誤差百分比

% 優化器參數設定 (Levenberg-Marquardt 演算法)
options = optimoptions('lsqnonlin', 'Display', 'none', 'Algorithm', 'levenberg-marquardt');

%% ============================================================
% 3. 進入熱力圖網格掃描
% ============================================================
total_points = numel(X_mesh);
fprintf('開始進行 TDOA 網格熱力圖掃描，總計點數：%d ...\n', total_points);

tic; % 開始計時
for i = 1:size(X_mesh, 1)
    for j = 1:size(X_mesh, 2)
        % 當前網格點的真實聲源位置
        true_x = X_mesh(i, j);
        true_y = Y_mesh(i, j);
        true_pos = [true_x; true_y];
        
        % ----------------------------------------------------
        % 正向模擬：計算理論 TOA 並引入 PRR 時間量化誤差
        % ----------------------------------------------------
        % 到各感測器的理論距離與到達時間 (TOA)
        distances = sqrt(sum((sensors - true_pos').^2, 2));
        trueArrivalTimes = distances / c;
        
        % 模擬硬體取樣率造成的量化誤差 (Quantization)
        quantizedArrivalTimes = round(trueArrivalTimes * PRR_Hz) / PRR_Hz;
        
        % 計算以 S1 為基準的時間差 TDOA 與距離差 d
        tdoa = quantizedArrivalTimes;
        tdoa(refIndex) = [];
        tdoa = tdoa - quantizedArrivalTimes(refIndex);
        d = c * tdoa;
        
        % ----------------------------------------------------
        % 逆向定位：使用 yn_TDOA 的核心架構求解
        % ----------------------------------------------------
        % 沒有 LLS 預估，每一次都固定從幾何質心 p0 = [0; 0] 出發
        [estimated_pos, ~] = lsqnonlin(@(p) tdoa_residuals_2d(p, sensors, d), p0, [], [], options);
        
        % ----------------------------------------------------
        % 誤差計算與百分比處理
        % ----------------------------------------------------
        % 計算 2D 幾何定位誤差 (公尺)
        error_m = sqrt(sum((estimated_pos - true_pos).^2));
        
        % 計算聲源到原點的真實距離 (用於計算百分比)
        true_dist_from_origin = sqrt(true_x^2 + true_y^2);
        
        if true_dist_from_origin < 1e-5
            error_percent = 0; % 原點處不計算除以零
        else
            error_percent = (error_m / true_dist_from_origin) * 100;
        end
        
        errorPercentMap(i, j) = error_percent;
    end
end
toc; % 結束計時

%% ============================================================
% 4. 誤差百分比分級 (全新指定的五大級距)
% ============================================================
classMap = zeros(size(errorPercentMap));
classMap(errorPercentMap < 1) = 1;                                 % 級別 1: <1%
classMap(errorPercentMap >= 1 & errorPercentMap < 5) = 2;          % 級別 2: 1 ~ 5%
classMap(errorPercentMap >= 5 & errorPercentMap < 10) = 3;         % 級別 3: 5 ~ 10%
classMap(errorPercentMap >= 10 & errorPercentMap < 20) = 4;        % 級別 4: 10 ~ 20%
classMap(errorPercentMap >= 20) = 5;                               % 級別 5: >20%

%% ============================================================
% 5. 熱力圖結果視覺化
% ============================================================
figure('Color', 'w', 'Position', [150, 120, 800, 650]);

% 繪製級距圖像
imageHandle = imagesc(xCenters, yCenters, classMap);
set(imageHandle, 'AlphaData', ~isnan(classMap));
set(gca, 'YDir', 'normal');
axis equal; axis tight; grid on; box on;

% 動態調整主要與次要網格線（維持畫面美觀與精準度）
ax = gca;
tick_step = half_side / 2;
ax.XTick = -half_side:tick_step:half_side;
ax.YTick = -half_side:tick_step:half_side;

ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.XAxis.MinorTickValues = -half_side:grid_spacing_m:half_side;
ax.YAxis.MinorTickValues = -half_side:grid_spacing_m:half_side;

% 套用經典 5 色 Colormap (由冷色到暖色，凸顯誤差嚴重性)
colormap([
    0.00, 0.45, 0.74;   % 級別 1: <1% (深藍 - 極準確)
    0.47, 0.67, 0.19;   % 級別 2: 1 ~ 5% (綠色 - 良好)
    0.93, 0.69, 0.13;   % 級別 3: 5 ~ 10% (黃色 - 輕微偏差)
    0.85, 0.33, 0.10;   % 級別 4: 10 ~ 20% (橘色 - 明顯誤差)
    0.49, 0.18, 0.56    % 級別 5: >20% (紫色 - 嚴重發散)
]);
caxis([0.5, 5.5]);

% 色標 TickLabels 渲染
cb = colorbar;
cb.Ticks = 1:5;
labelDash = char(8211);
cb.TickLabels = { ...
    '<5%', ...
    ['5' labelDash '10%'], ...
    ['10' labelDash '20%'], ...
    ['20' labelDash '100%'], ...
    '>100%' ...
};
ylabel(cb, '定位誤差百分比 (Error %)', 'FontSize', 11, 'FontWeight', 'bold');

hold on;

% 繪製 5 顆感測器位置
plot(sensors(:,1), sensors(:,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 8, 'LineWidth', 2);

% 圖表標題與標籤
xlabel('X 軸位置 (公尺)', 'FontSize', 11);
ylabel('Y 軸位置 (公尺)', 'FontSize', 11);
title_str = sprintf('2D TDOA 誤差百分比熱力圖\n 範圍: %d m x %d m  |  PRR: %d Hz', area_side_m, area_side_m, PRR_Hz);
title(title_str, 'FontSize', 12, 'FontWeight', 'bold');

hold off;

%% ============================================================
% 6. 輔助子函數：2D TDOA 殘差計算 (Local Function)
% ============================================================
function res = tdoa_residuals_2d(p, sensors, d)
    % p: [x; y] 估計的聲源位置 (米)
    % sensors: Nx2 感測器位置矩陣 (米)
    % d: (N-1)x1 到達距離差矩陣 (米)
    x = p(1); y = p(2);
    
    % 到參考感測器 S1 的距離
    r0 = sqrt((x - sensors(1,1))^2 + (y - sensors(1,2))^2);
    res = zeros(length(d), 1);
    
    % 計算其餘感測器相對於 S1 的殘差
    for k = 2:size(sensors, 1)
        rk = sqrt((x - sensors(k,1))^2 + (y - sensors(k,2))^2);
        res(k-1) = rk - r0 - d(k-1); 
    end
end