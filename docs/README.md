# TDOA Source Localization Method

本文件整理 TDOA（Time Difference of Arrival）聲源定位方法。  
定位流程分成兩階段：

1. 使用線性 Least Squares 求初始位置。
2. 使用 Levenberg-Marquardt Method 對非線性 TDOA 方程式做迭代修正。

---

## 1. Problem Definition

未知聲源位置為：

$$
p=
\begin{bmatrix}
x\\
y
\end{bmatrix}
$$

第 $i$ 顆感測器位置為：

$$
s_i=
\begin{bmatrix}
x_i\\
y_i
\end{bmatrix}
$$

reference sensor 位置為：

$$
s_{ref}=
\begin{bmatrix}
x_{ref}\\
y_{ref}
\end{bmatrix}
$$

聲源到第 $i$ 顆感測器的距離：

$$
d_i=\|p-s_i\|=\sqrt{(x-x_i)^2+(y-y_i)^2}
$$

聲源到 reference sensor 的距離：

$$
d_{ref}=\|p-s_{ref}\|=\sqrt{(x-x_{ref})^2+(y-y_{ref})^2}
$$

---

## 2. TDOA Measurement Model

第 $i$ 顆感測器接收到訊號的時間為 $t_i$，reference sensor 接收到訊號的時間為 $t_{ref}$。

定義時間差：

$$
\Delta t_i=t_i-t_{ref}
$$

乘上聲速 $c$，得到距離差：

$$
\rho_i=c\Delta t_i=c(t_i-t_{ref})
$$

本文採用的 TDOA 定義為：

$$
\rho_i=d_i-d_{ref}
$$

因此核心方程式為：

$$
\|p-s_i\|-\|p-s_{ref}\|=\rho_i
$$

展開後：

$$
\sqrt{(x-x_i)^2+(y-y_i)^2}-\sqrt{(x-x_{ref})^2+(y-y_{ref})^2}=\rho_i
$$

這是 TDOA 定位的原始非線性方程式。

---

## 3. Note on Signal Delay Sign

若聲源訊號為：

$$
s(t)=A\cos(2\pi ft)
$$

延遲 $\tau_i$ 後到達第 $i$ 顆感測器，通常寫成：

$$
s_i(t)=A\cos(2\pi f(t-\tau_i))
$$

---

## 4. Issue in the Original Linear Formula

原投影片將方程式整理成：

$$
AX=b
$$

其中：

$$
X=
\begin{bmatrix}
x\\
y
\end{bmatrix}
$$

且：

$$
b_i=
x_{ref}^2-x_i^2+y_{ref}^2-y_i^2+d_i^2-d_{ref}^2
$$

這裡的問題是：TDOA 量測只能得到距離差

$$
\rho_i=d_i-d_{ref}
$$

但無法直接得到 $d_i$ 和 $d_{ref}$ 的個別數值。

因此：

$$
d_i^2-d_{ref}^2
$$

不是已知量，不能直接放進 $b_i$。

所以線性 Least Squares 不能只解：

$$
X=
\begin{bmatrix}
x\\
y
\end{bmatrix}
$$

必須額外引入未知數：

$$
r=d_{ref}
$$

因此正確的未知向量應為：

$$
u=
\begin{bmatrix}
x\\
y\\
r
\end{bmatrix}
$$

---

## 5. Correct Linear Least Squares Form

由 TDOA 定義：

$$
\rho_i=d_i-d_{ref}
$$

令：

$$
r=d_{ref}
$$

則：

$$
d_i=r+\rho_i
$$

兩邊平方：

$$
d_i^2=(r+\rho_i)^2
$$

展開：

$$
d_i^2=r^2+2\rho_i r+\rho_i^2
$$

由於 $r=d_{ref}$，所以：

$$
r^2=d_{ref}^2
$$

因此：

$$
d_i^2=d_{ref}^2+2\rho_i r+\rho_i^2
$$

移項：

$$
d_i^2-d_{ref}^2=2\rho_i r+\rho_i^2
$$

---

## 6. Distance-Square Expansion

距離平方為：

$$
d_i^2=(x-x_i)^2+(y-y_i)^2
$$

$$
d_{ref}^2=(x-x_{ref})^2+(y-y_{ref})^2
$$

因此：

$$
d_i^2-d_{ref}^2=(x-x_i)^2+(y-y_i)^2-(x-x_{ref})^2-(y-y_{ref})^2
$$

展開：

$$
d_i^2-d_{ref}^2=2x(x_{ref}-x_i)+2y(y_{ref}-y_i)+x_i^2+y_i^2-x_{ref}^2-y_{ref}^2
$$

又因為：

$$
d_i^2-d_{ref}^2=2\rho_i r+\rho_i^2
$$

所以：

$$
2x(x_{ref}-x_i)+2y(y_{ref}-y_i)+x_i^2+y_i^2-x_{ref}^2-y_{ref}^2=2\rho_i r+\rho_i^2
$$

整理成線性形式：

$$
2(x_{ref}-x_i)x+2(y_{ref}-y_i)y-2\rho_i r=x_{ref}^2+y_{ref}^2-x_i^2-y_i^2+\rho_i^2
$$

---

## 7. Matrix Form

定義：

$$
u=
\begin{bmatrix}
x\\
y\\
r
\end{bmatrix}
$$

則對第 $i$ 顆感測器：

$$
A_i=
\begin{bmatrix}
2(x_{ref}-x_i) & 2(y_{ref}-y_i) & -2\rho_i
\end{bmatrix}
$$

$$
b_i=
x_{ref}^2+y_{ref}^2-x_i^2-y_i^2+\rho_i^2
$$

整體系統為：

$$
Au=b
$$

其中：

$$
A=
\begin{bmatrix}
2(x_{ref}-x_1) & 2(y_{ref}-y_1) & -2\rho_1\\
2(x_{ref}-x_2) & 2(y_{ref}-y_2) & -2\rho_2\\
\vdots & \vdots & \vdots\\
2(x_{ref}-x_n) & 2(y_{ref}-y_n) & -2\rho_n
\end{bmatrix}
$$

$$
b=
\begin{bmatrix}
x_{ref}^2+y_{ref}^2-x_1^2-y_1^2+\rho_1^2\\
x_{ref}^2+y_{ref}^2-x_2^2-y_2^2+\rho_2^2\\
\vdots\\
x_{ref}^2+y_{ref}^2-x_n^2-y_n^2+\rho_n^2
\end{bmatrix}
$$

---

## 8. Least Squares Solution

當感測器數量大於未知數數量時，系統通常為 overdetermined system：

$$
Au\approx b
$$

使用 Least Squares：

$$
\min_u \|Au-b\|^2
$$

normal equation 為：

$$
A^TAu=A^Tb
$$

理論解：

$$
u=(A^TA)^{-1}A^Tb
$$


MATLAB：

```matlab
u = A \ b;
```

Python：

```python
u = np.linalg.lstsq(A, b, rcond=None)[0]
```

取得初始位置：

$$
p_0=\begin{bmatrix}u_1\\u_2\end{bmatrix}=\begin{bmatrix}x\\y\end{bmatrix}
$$

---

## 9. Nonlinear Residual Function

Least Squares 只用來取得初始值。  
接著回到原始 TDOA 非線性方程式進行修正。

本文定義：

$$
\rho_i=d_i-d_{ref}
$$

所以第 $i$ 顆感測器的殘差為：

$$
e_i(p)=\|p-s_i\|-\|p-s_{ref}\|-\rho_i
$$

展開：

$$
e_i(p)=\sqrt{(x-x_i)^2+(y-y_i)^2}-\sqrt{(x-x_{ref})^2+(y-y_{ref})^2}-\rho_i
$$

所有殘差組成：

$$
e(p)=
\begin{bmatrix}
e_1(p)\\
e_2(p)\\
\vdots\\
e_n(p)
\end{bmatrix}
$$

非線性最佳化目標：

$$
\min_p \frac{1}{2}\|e(p)\|^2
$$

---

## 10. Jacobian Matrix

對殘差函數：

$$
e_i(p)=\|p-s_i\|-\|p-s_{ref}\|-\rho_i
$$

計算偏微分。

對 $x$：

$$
\frac{\partial e_i}{\partial x}=\frac{x-x_i}{\|p-s_i\|}-\frac{x-x_{ref}}{\|p-s_{ref}\|}
$$

對 $y$：

$$
\frac{\partial e_i}{\partial y}=\frac{y-y_i}{\|p-s_i\|}-\frac{y-y_{ref}}{\|p-s_{ref}\|}
$$

因此：

$$
J_i(p)=\begin{bmatrix}\frac{x-x_i}{\|p-s_i\|}-\frac{x-x_{ref}}{\|p-s_{ref}\|}&\frac{y-y_i}{\|p-s_i\|}-\frac{y-y_{ref}}{\|p-s_{ref}\|}\end{bmatrix}
$$

---

## 11. Levenberg-Marquardt Update

在第 $k$ 次迭代，位置為：

$$
p_k
$$

更新量為：

$$
\Delta p_k
$$

一階近似：

$$
e(p_k+\Delta p_k)
\approx
e(p_k)+J(p_k)\Delta p_k
$$

令：

$$
e_k=e(p_k)
$$

$$
J_k=J(p_k)
$$

LM 要解：

$$
\min_{\Delta p_k}
\frac{1}{2}\|e_k+J_k\Delta p_k\|^2
+
\frac{\lambda}{2}\|\Delta p_k\|^2
$$

對 $\Delta p_k$ 微分並令為 0：

$$
(J_k^TJ_k+\lambda I)\Delta p_k=-J_k^Te_k
$$

因此：

$$
\Delta p_k=-(J_k^TJ_k+\lambda I)^{-1}J_k^Te_k
$$

實作時建議解線性系統。

MATLAB：

```matlab
delta_p = -(J.' * J + lambda * eye(size(J,2))) \ (J.' * e);
```

Python：

```python
delta_p = -np.linalg.solve(J.T @ J + lambda * np.eye(J.shape[1]), J.T @ e)
```

更新位置：

$$
p_{k+1}=p_k+\Delta p_k
$$

---

## 12. 3D Extension

若改為三維定位：

$$
p=
\begin{bmatrix}
x\\
y\\
z
\end{bmatrix}
$$

距離為：

$$
d_i=\sqrt{(x-x_i)^2+(y-y_i)^2+(z-z_i)^2}
$$

TDOA 殘差仍為：

$$
e_i(p)=\|p-s_i\|-\|p-s_{ref}\|-\rho_i
$$

Jacobian 變為：

$$
J_i(p)=\begin{bmatrix}\frac{x-x_i}{\|p-s_i\|}-\frac{x-x_{ref}}{\|p-s_{ref}\|}&\frac{y-y_i}{\|p-s_i\|}-\frac{y-y_{ref}}{\|p-s_{ref}\|}&\frac{z-z_i}{\|p-s_i\|}-\frac{z-z_{ref}}{\|p-s_{ref}\|}\end{bmatrix}
$$

三維線性初始解的未知向量為：

$$
u=
\begin{bmatrix}
x\\
y\\
z\\
r
\end{bmatrix}
$$

其中：

$$
r=d_{ref}
$$

第 $i$ 列矩陣為：

$$
A_i=\begin{bmatrix}2(x_{ref}-x_i) &2(y_{ref}-y_i) &2(z_{ref}-z_i) &-2\rho_i\end{bmatrix}
$$

$$
b_i=x_{ref}^2+y_{ref}^2+z_{ref}^2-x_i^2-y_i^2-z_i^2+\rho_i^2
$$

---

## 13. Algorithm Flow

1. 設定感測器座標 $s_1,s_2,\dots,s_N$。
2. 選擇 reference sensor：$s_{ref}$。
3. 量測時間差：

$$
\Delta t_i=t_i-t_{ref}
$$

4. 轉為距離差：

$$
\rho_i=c\Delta t_i
$$

5. 建立線性系統：

$$
Au=b
$$

6. 求 Least Squares 初始解：

$$
u=A\backslash b
$$

7. 取初始位置：

$$
p_0=
\begin{bmatrix}
x\\
y
\end{bmatrix}
$$

8. 建立非線性殘差：

$$
e_i(p)=\|p-s_i\|-\|p-s_{ref}\|-\rho_i
$$

9. 計算 Jacobian $J(p)$。

10. 使用 LM 更新：

$$
(J_k^TJ_k+\lambda I)\Delta p_k=-J_k^Te_k
$$

11. 更新位置：

$$
p_{k+1}=p_k+\Delta p_k
$$

12. 當 $\|\Delta p_k\|$ 足夠小或殘差不再下降時停止。

---

## 14. Sign Convention

本文使用：

$$
\rho_i=c(t_i-t_{ref})=d_i-d_{ref}
$$

因此殘差必須寫成：

$$
e_i(p)=\|p-s_i\|-\|p-s_{ref}\|-\rho_i
$$

若改用：

$$
\rho_i=c(t_{ref}-t_i)=d_{ref}-d_i
$$

則殘差要同步改成：

$$
e_i(p)=\|p-s_{ref}\|-\|p-s_i\|-\rho_i
$$

TDOA 實作時最重要的是正負號必須全程一致。

---

## 15. Summary

TDOA 核心量測式：

$$
\rho_i=c(t_i-t_{ref})=d_i-d_{ref}
$$

核心非線性模型：

$$
\|p-s_i\|-\|p-s_{ref}\|=\rho_i
$$

線性初始解：

$$
Au=b
$$

其中二維未知量為：

$$
u=
\begin{bmatrix}
x\\
y\\
r
\end{bmatrix}
$$

LM 非線性修正：

$$
(J_k^TJ_k+\lambda I)\Delta p_k=-J_k^Te_k
$$

最終更新：

$$
p_{k+1}=p_k+\Delta p_k
$$
