# TDOA Source Localization Method

本文件整理以到達時間差（Time Difference of Arrival, TDOA）進行聲源定位的基本模型與求解方法。在已知感測器位置及訊號傳播速度的條件下，定位流程可分為兩個階段：

1. 使用線性最小平方法（Least Squares, LS）估計聲源的初始位置。
2. 使用 Levenberg–Marquardt（LM）演算法求解非線性 TDOA 方程式，進一步修正定位結果。

以下先說明二維定位模型，並於後文延伸至三維情況。

---

## 1. 問題定義（Problem Definition）

令未知聲源位置為

$$
p=
\begin{bmatrix}
x\\
y
\end{bmatrix},
$$

第 $i$ 顆感測器的位置為

$$
s_i=
\begin{bmatrix}
x_i\\
y_i
\end{bmatrix},
$$

參考感測器（reference sensor）的位置為

$$
s_{\mathrm{ref}}=
\begin{bmatrix}
x_{\mathrm{ref}}\\
y_{\mathrm{ref}}
\end{bmatrix}.
$$

聲源至第 $i$ 顆感測器的距離為

$$
d_i=\lVert p-s_i\rVert_2
=\sqrt{(x-x_i)^2+(y-y_i)^2},
$$

聲源至參考感測器的距離為

$$
d_{\mathrm{ref}}=\lVert p-s_{\mathrm{ref}}\rVert_2
=\sqrt{(x-x_{\mathrm{ref}})^2+(y-y_{\mathrm{ref}})^2}.
$$

---

## 2. TDOA 量測模型（Measurement Model）

令第 $i$ 顆感測器與參考感測器接收到訊號的時間分別為 $t_i$ 與 $t_{\mathrm{ref}}$。兩者的到達時間差定義為

$$
\Delta t_i=t_i-t_{\mathrm{ref}}.
$$

將時間差乘以訊號傳播速度 $c$，可得到距離差

$$
\rho_i=c\Delta t_i=c(t_i-t_{\mathrm{ref}}).
$$

本文件採用下列正負號定義：

$$
\rho_i=d_i-d_{\mathrm{ref}}.
$$

因此，TDOA 定位的核心方程式為

$$
\lVert p-s_i\rVert_2-\lVert p-s_{\mathrm{ref}}\rVert_2=\rho_i.
$$

展開後可寫成

$$
\sqrt{(x-x_i)^2+(y-y_i)^2}
-\sqrt{(x-x_{\mathrm{ref}})^2+(y-y_{\mathrm{ref}})^2}
=\rho_i.
$$

由於未知位置 $(x,y)$ 位於平方根內，因此這是一組非線性方程式。

---

## 3. 訊號延遲與正負號（Signal Delay Sign）

若原始訊號為

$$
s(t)=A\cos(2\pi ft),
$$

訊號延遲 $\tau_i$ 後到達第 $i$ 顆感測器，可表示為

$$
s_i(t)=A\cos\left[2\pi f(t-\tau_i)\right].
$$

因此，較晚接收到訊號的感測器具有較大的到達時間。實作時必須讓時間差、距離差與殘差函數採用相同的正負號定義。

---

## 4. 常見的線性化問題

若直接將定位方程式整理為

$$
AX=b,
$$

並令

$$
X=
\begin{bmatrix}
x\\
y
\end{bmatrix},
$$

可能會在 $b_i$ 中出現

$$
d_i^2-d_{\mathrm{ref}}^2.
$$

然而，TDOA 量測只能取得距離差

$$
\rho_i=d_i-d_{\mathrm{ref}},
$$

無法分別得知 $d_i$ 與 $d_{\mathrm{ref}}$。因此，$d_i^2-d_{\mathrm{ref}}^2$ 不能直接視為已知量。

為了建立可求解的線性系統，需引入輔助未知數

$$
r=d_{\mathrm{ref}},
$$

並將未知向量改寫為

$$
u=
\begin{bmatrix}
x\\
y\\
r
\end{bmatrix}.
$$

---

## 5. 線性最小平方法的推導

由

$$
\rho_i=d_i-d_{\mathrm{ref}}
$$

及 $r=d_{\mathrm{ref}}$，可得

$$
d_i=r+\rho_i.
$$

兩邊平方後得到

$$
d_i^2=(r+\rho_i)^2
=r^2+2\rho_i r+\rho_i^2.
$$

由於 $r^2=d_{\mathrm{ref}}^2$，因此

$$
d_i^2-d_{\mathrm{ref}}^2=2\rho_i r+\rho_i^2.
$$

另一方面，距離平方差可展開為

$$
\begin{aligned}
d_i^2-d_{\mathrm{ref}}^2
={}&(x-x_i)^2+(y-y_i)^2\\
&-(x-x_{\mathrm{ref}})^2-(y-y_{\mathrm{ref}})^2.
\end{aligned}
$$

整理後可得

$$
\begin{aligned}
d_i^2-d_{\mathrm{ref}}^2
={}&2x(x_{\mathrm{ref}}-x_i)
+2y(y_{\mathrm{ref}}-y_i)\\
&+x_i^2+y_i^2-x_{\mathrm{ref}}^2-y_{\mathrm{ref}}^2.
\end{aligned}
$$

將兩式聯立，最終可整理為線性形式

$$
\begin{aligned}
2(x_{\mathrm{ref}}-x_i)x
&+2(y_{\mathrm{ref}}-y_i)y
-2\rho_i r\\
&=x_{\mathrm{ref}}^2+y_{\mathrm{ref}}^2-x_i^2-y_i^2+\rho_i^2.
\end{aligned}
$$

---

## 6. 矩陣形式（Matrix Form）

對每一顆非參考感測器，定義

$$
A_i=
\begin{bmatrix}
2(x_{\mathrm{ref}}-x_i) &
2(y_{\mathrm{ref}}-y_i) &
-2\rho_i
\end{bmatrix},
$$

以及

$$
b_i=x_{\mathrm{ref}}^2+y_{\mathrm{ref}}^2-x_i^2-y_i^2+\rho_i^2.
$$

將所有非參考感測器的方程式堆疊後，可得

$$
Au=b,
$$

其中

$$
u=
\begin{bmatrix}
x\\
y\\
r
\end{bmatrix}.
$$

若共有 $M$ 組相對於參考感測器的 TDOA 量測，則 $A\in\mathbb{R}^{M\times3}$、$u\in\mathbb{R}^{3}$、$b\in\mathbb{R}^{M}$。

---

## 7. 最小平方法求初始位置（Least Squares Solution）

當獨立量測方程式多於未知數時，系統通常為超定系統（overdetermined system）：

$$
Au\approx b.
$$

此時可透過最小化平方殘差求得線性解：

$$
\min_u \lVert Au-b\rVert_2^2.
$$

其常態方程式（normal equation）為

$$
A^TAu=A^Tb.
$$

若 $A^TA$ 可逆，理論解為

$$
u=(A^TA)^{-1}A^Tb.
$$

實際計算時不建議直接求反矩陣，應使用數值較穩定的線性求解器或虛擬反矩陣。

MATLAB：

```matlab
u = A \ b;
```

Python：

```python
u = np.linalg.lstsq(A, b, rcond=None)[0]
```

由線性解取得 LM 演算法的初始位置：

$$
p_0=
\begin{bmatrix}
u_1\\
u_2
\end{bmatrix}
=
\begin{bmatrix}
x\\
y
\end{bmatrix}.
$$

在二維定位中，未知向量 $u=[x,y,r]^T$ 共有三個分量，因此至少需要三組獨立的 TDOA 方程式，亦即至少四顆感測器（包含參考感測器）。若方程式不足或矩陣秩不足，則無法由此線性系統取得唯一初始解，需另行設定 $p_0$。

---

## 8. 非線性殘差函數（Nonlinear Residual Function）

線性最小平方法主要用於取得初始位置。接著回到原始 TDOA 模型，使用非線性最佳化進一步修正定位結果。

依照本文件的正負號定義，第 $i$ 組量測的殘差為

$$
e_i(p)=\lVert p-s_i\rVert_2
-\lVert p-s_{\mathrm{ref}}\rVert_2
-\rho_i.
$$

展開後為

$$
\begin{aligned}
e_i(p)={}&\sqrt{(x-x_i)^2+(y-y_i)^2}\\
&-\sqrt{(x-x_{\mathrm{ref}})^2+(y-y_{\mathrm{ref}})^2}
-\rho_i.
\end{aligned}
$$

將所有量測殘差組成向量

$$
e(p)=
\begin{bmatrix}
e_1(p)\\
e_2(p)\\
\vdots\\
e_M(p)
\end{bmatrix},
$$

非線性最佳化的目標為

$$
\min_p \frac{1}{2}\lVert e(p)\rVert_2^2.
$$

---

## 9. Jacobian 矩陣

對殘差函數

$$
e_i(p)=\lVert p-s_i\rVert_2
-\lVert p-s_{\mathrm{ref}}\rVert_2
-\rho_i,
$$

分別對 $x$ 與 $y$ 計算偏微分：

$$
\frac{\partial e_i}{\partial x}
=\frac{x-x_i}{\lVert p-s_i\rVert_2}
-\frac{x-x_{\mathrm{ref}}}{\lVert p-s_{\mathrm{ref}}\rVert_2},
$$

$$
\frac{\partial e_i}{\partial y}
=\frac{y-y_i}{\lVert p-s_i\rVert_2}
-\frac{y-y_{\mathrm{ref}}}{\lVert p-s_{\mathrm{ref}}\rVert_2}.
$$

因此，第 $i$ 列 Jacobian 可寫為

$$
J_i(p)=
\begin{bmatrix}
\dfrac{x-x_i}{\lVert p-s_i\rVert_2}
-\dfrac{x-x_{\mathrm{ref}}}{\lVert p-s_{\mathrm{ref}}\rVert_2}
&
\dfrac{y-y_i}{\lVert p-s_i\rVert_2}
-\dfrac{y-y_{\mathrm{ref}}}{\lVert p-s_{\mathrm{ref}}\rVert_2}
\end{bmatrix}.
$$

---

## 10. Levenberg–Marquardt 更新

在第 $k$ 次迭代中，令目前的位置估計為 $p_k$，更新量為 $\Delta p_k$。對殘差函數進行一階近似：

$$
e(p_k+\Delta p_k)
\approx e(p_k)+J(p_k)\Delta p_k.
$$

令

$$
e_k=e(p_k),\qquad J_k=J(p_k),
$$

則 LM 演算法求解下列阻尼最小平問題：

$$
\min_{\Delta p_k}
\left(
\frac{1}{2}\lVert e_k+J_k\Delta p_k\rVert_2^2
+\frac{\lambda}{2}\lVert\Delta p_k\rVert_2^2
\right).
$$

對 $\Delta p_k$ 微分並令其為零，可得

$$
(J_k^TJ_k+\lambda I)\Delta p_k=-J_k^Te_k.
$$

因此，更新量可表示為

$$
\Delta p_k
=-(J_k^TJ_k+\lambda I)^{-1}J_k^Te_k.
$$

實作時同樣應直接求解線性系統，而非顯式計算反矩陣。

MATLAB：

```matlab
delta_p = -(J.' * J + lambda * eye(size(J, 2))) \ (J.' * e);
```

Python：

```python
delta_p = -np.linalg.solve(
    J.T @ J + lambda * np.eye(J.shape[1]),
    J.T @ e,
)
```

最後更新聲源位置：

$$
p_{k+1}=p_k+\Delta p_k.
$$

---

## 11. 三維定位延伸（3D Extension）

若將模型延伸至三維空間，聲源位置改寫為

$$
p=
\begin{bmatrix}
x\\
y\\
z
\end{bmatrix},
$$

聲源至第 $i$ 顆感測器的距離為

$$
d_i=\sqrt{(x-x_i)^2+(y-y_i)^2+(z-z_i)^2}.
$$

TDOA 殘差形式維持不變：

$$
e_i(p)=\lVert p-s_i\rVert_2
-\lVert p-s_{\mathrm{ref}}\rVert_2
-\rho_i.
$$

三維 Jacobian 的第 $i$ 列為

$$
J_i(p)=
\begin{bmatrix}
\dfrac{x-x_i}{\lVert p-s_i\rVert_2}
-\dfrac{x-x_{\mathrm{ref}}}{\lVert p-s_{\mathrm{ref}}\rVert_2}
&
\dfrac{y-y_i}{\lVert p-s_i\rVert_2}
-\dfrac{y-y_{\mathrm{ref}}}{\lVert p-s_{\mathrm{ref}}\rVert_2}
&
\dfrac{z-z_i}{\lVert p-s_i\rVert_2}
-\dfrac{z-z_{\mathrm{ref}}}{\lVert p-s_{\mathrm{ref}}\rVert_2}
\end{bmatrix}.
$$

三維線性初始解的未知向量為

$$
u=
\begin{bmatrix}
x\\
y\\
z\\
r
\end{bmatrix},
\qquad r=d_{\mathrm{ref}}.
$$

對應的第 $i$ 列矩陣與右側向量分量為

$$
A_i=
\begin{bmatrix}
2(x_{\mathrm{ref}}-x_i) &
2(y_{\mathrm{ref}}-y_i) &
2(z_{\mathrm{ref}}-z_i) &
-2\rho_i
\end{bmatrix},
$$

$$
b_i=x_{\mathrm{ref}}^2+y_{\mathrm{ref}}^2+z_{\mathrm{ref}}^2
-x_i^2-y_i^2-z_i^2+\rho_i^2.
$$

---

## 12. 演算法流程（Algorithm Flow）

1. 設定所有感測器的位置 $s_1,s_2,\ldots,s_N$。
2. 選定參考感測器 $s_{\mathrm{ref}}$。
3. 計算各感測器相對於參考感測器的到達時間差：

   $$
   \Delta t_i=t_i-t_{\mathrm{ref}}.
   $$

4. 將時間差轉換為距離差：

   $$
   \rho_i=c\Delta t_i.
   $$

5. 建立線性系統 $Au=b$，並以最小平方法求得初始位置 $p_0$。
6. 根據原始 TDOA 模型建立非線性殘差 $e(p)$ 與 Jacobian $J(p)$。
7. 使用 LM 演算法求得更新量：

   $$
   (J_k^TJ_k+\lambda I)\Delta p_k=-J_k^Te_k.
   $$

8. 更新位置估計：

   $$
   p_{k+1}=p_k+\Delta p_k.
   $$

9. 當 $\lVert\Delta p_k\rVert_2$ 足夠小、殘差不再明顯下降，或達到最大迭代次數時停止。

---

## 13. 正負號定義（Sign Convention）

本文件採用

$$
\rho_i=c(t_i-t_{\mathrm{ref}})=d_i-d_{\mathrm{ref}},
$$

因此殘差必須寫成

$$
e_i(p)=\lVert p-s_i\rVert_2
-\lVert p-s_{\mathrm{ref}}\rVert_2
-\rho_i.
$$

若改用

$$
\rho_i=c(t_{\mathrm{ref}}-t_i)=d_{\mathrm{ref}}-d_i,
$$

則殘差也必須同步改為

$$
e_i(p)=\lVert p-s_{\mathrm{ref}}\rVert_2
-\lVert p-s_i\rVert_2
-\rho_i.
$$

TDOA 實作的關鍵並非採用哪一種定義，而是時間差、距離差及殘差函數的正負號必須全程一致。

---

## 14. 總結（Summary）

TDOA 的距離差量測式為

$$
\rho_i=c(t_i-t_{\mathrm{ref}})=d_i-d_{\mathrm{ref}}.
$$

核心非線性定位模型為

$$
\lVert p-s_i\rVert_2
-\lVert p-s_{\mathrm{ref}}\rVert_2
=\rho_i.
$$

引入輔助未知數 $r=d_{\mathrm{ref}}$ 後，可建立線性系統

$$
Au=b
$$

以取得初始位置，再使用 LM 演算法進行非線性修正：

$$
(J_k^TJ_k+\lambda I)\Delta p_k=-J_k^Te_k,
$$

$$
p_{k+1}=p_k+\Delta p_k.
$$

此兩階段方法結合線性初始估計與非線性迭代修正，可降低初始值選擇對最終定位結果的影響。
