# TDOA 聲源定位

使用 MATLAB，根據訊號到達時間差估計二維聲源位置。

| 檔案／資料夾 | 內容 |
| --- | --- |
| [TDOA_math.md](TDOA_math.md) | TDOA 數學推導（目前為 3D 筆記）。 |
| [tdoa_localization_2d.m](tdoa_localization_2d.m) | 2D 定位演算法，使用 Linear LS 初始解與 LM 修正；三顆感測器直接以原點初始化。 |
| [data/](data/) | TOA 測資，待上傳。 |
| [results/](results/) | 定位結果圖，待上傳。 |

執行程式需要 MATLAB 與 Optimization Toolbox；目前測資設定於程式內。

