# 3D TDOA Source Localization Notes

## 1. 基本設定

未知聲源位置：

```math
p=
\begin{bmatrix}
x\\
y\\
z
\end{bmatrix}
```

第 $`i`$ 顆感測器位置：

```math
s_i=
\begin{bmatrix}
x_i\\
y_i\\
z_i
\end{bmatrix}
```

參考感測器位置：

```math
s_{\mathrm{ref}}=
\begin{bmatrix}
x_{\mathrm{ref}}\\
y_{\mathrm{ref}}\\
z_{\mathrm{ref}}
\end{bmatrix}
```

聲源至第 $`i`$ 顆感測器的距離：

```math
d_i
=
\lVert p-s_i\rVert_2
=
\sqrt{
(x-x_i)^2+
(y-y_i)^2+
(z-z_i)^2
}
```

聲源至參考感測器的距離：

```math
d_{\mathrm{ref}}
=
\lVert p-s_{\mathrm{ref}}\rVert_2
=
\sqrt{
(x-x_{\mathrm{ref}})^2+
(y-y_{\mathrm{ref}})^2+
(z-z_{\mathrm{ref}})^2
}
```

---

## 2. TDOA 量測模型

第 $`i`$ 顆感測器與參考感測器的訊號到達時間分別為 $`t_i`$ 與 $`t_{\mathrm{ref}}`$。

到達時間差：

```math
\Delta t_i=t_i-t_{\mathrm{ref}}
```

將時間差乘以訊號傳播速度 $`c`$，得到距離差：

```math
\rho_i
=
c\Delta t_i
=
c(t_i-t_{\mathrm{ref}})
```

對應的距離關係：

```math
\rho_i=d_i-d_{\mathrm{ref}}
```

TDOA 定位方程式：

```math
\lVert p-s_i\rVert_2
-
\lVert p-s_{\mathrm{ref}}\rVert_2
=
\rho_i
```

展開後：

```math
\sqrt{
(x-x_i)^2+
(y-y_i)^2+
(z-z_i)^2
}
-
\sqrt{
(x-x_{\mathrm{ref}})^2+
(y-y_{\mathrm{ref}})^2+
(z-z_{\mathrm{ref}})^2
}
=
\rho_i
```

未知位置 $`(x,y,z)`$ 位於平方根內，因此需要使用非線性方法求解。

---

## 3. 線性初始解

TDOA 只能取得距離差 $`\rho_i`$，無法分別取得 $`d_i`$ 與 $`d_{\mathrm{ref}}`$。

引入輔助未知數：

```math
r=d_{\mathrm{ref}}
```

未知向量改寫為：

```math
u=
\begin{bmatrix}
x\\
y\\
z\\
r
\end{bmatrix}
```

由

```math
\rho_i=d_i-d_{\mathrm{ref}}
```

可得

```math
d_i=r+\rho_i
```

兩邊平方：

```math
d_i^2
=
(r+\rho_i)^2
=
r^2+2\rho_i r+\rho_i^2
```

因為

```math
r^2=d_{\mathrm{ref}}^2
```

所以

```math
d_i^2-d_{\mathrm{ref}}^2
=
2\rho_i r+\rho_i^2
```

距離平方差也可寫成：

```math
\begin{aligned}
d_i^2-d_{\mathrm{ref}}^2
={}&
(x-x_i)^2+
(y-y_i)^2+
(z-z_i)^2\\
&-
(x-x_{\mathrm{ref}})^2
-
(y-y_{\mathrm{ref}})^2
-
(z-z_{\mathrm{ref}})^2
\end{aligned}
```

展開並整理：

```math
\begin{aligned}
d_i^2-d_{\mathrm{ref}}^2
={}&
2x(x_{\mathrm{ref}}-x_i)
+
2y(y_{\mathrm{ref}}-y_i)\\
&+
2z(z_{\mathrm{ref}}-z_i)
+
x_i^2+y_i^2+z_i^2\\
&-
x_{\mathrm{ref}}^2
-
y_{\mathrm{ref}}^2
-
z_{\mathrm{ref}}^2
\end{aligned}
```

與 $`2\rho_i r+\rho_i^2`$ 聯立後：

```math
\begin{aligned}
2(x_{\mathrm{ref}}-x_i)x
&+
2(y_{\mathrm{ref}}-y_i)y
+
2(z_{\mathrm{ref}}-z_i)z\\
&-
2\rho_i r\\
={}&
x_{\mathrm{ref}}^2
+
y_{\mathrm{ref}}^2
+
z_{\mathrm{ref}}^2\\
&-
x_i^2
-
y_i^2
-
z_i^2
+
\rho_i^2
\end{aligned}
```

---

## 4. 矩陣形式

對每一顆非參考感測器，定義：

```math
A_i=
\begin{bmatrix}
2(x_{\mathrm{ref}}-x_i) &
2(y_{\mathrm{ref}}-y_i) &
2(z_{\mathrm{ref}}-z_i) &
-2\rho_i
\end{bmatrix}
```

```math
b_i=
x_{\mathrm{ref}}^2
+
y_{\mathrm{ref}}^2
+
z_{\mathrm{ref}}^2
-
x_i^2
-
y_i^2
-
z_i^2
+
\rho_i^2
```

將所有非參考感測器的方程式排列後：

```math
Au=b
```

若共有 $`N`$ 顆感測器，其中一顆為參考感測器，則：

```math
A\in\mathbb{R}^{(N-1)\times4}
```

```math
u\in\mathbb{R}^{4}
```

```math
b\in\mathbb{R}^{N-1}
```

三維線性定位共有四個未知量 $`x`$、$`y`$、$`z`$ 與 $`r`$，因此至少需要四組獨立的 TDOA 方程式，也就是至少五顆感測器。

另外需滿足：

```math
\operatorname{rank}(A)=4
```

若感測器配置使矩陣秩不足，便無法得到唯一的三維線性初始解。

---

## 5. Least Squares 初始位置

當 TDOA 方程式數量多於未知數數量時，可使用最小平方法求解：

```math
\min_u \lVert Au-b\rVert_2^2
```

對應的常態方程式：

```math
A^TAu=A^Tb
```

理論解：

```math
u=(A^TA)^{-1}A^Tb
```

實際計算時，不直接計算反矩陣，而是使用線性求解器或虛擬反矩陣。

MATLAB：

```matlab
u = A \ b;
```

Python：

```python
u = np.linalg.lstsq(A, b, rcond=None)[0]
```

取出位置分量，作為非線性最佳化的初始位置：

```math
p_0=
\begin{bmatrix}
u_1\\
u_2\\
u_3
\end{bmatrix}
=
\begin{bmatrix}
x\\
y\\
z
\end{bmatrix}
```

---

## 6. 非線性殘差函數

線性最小平方法得到初始位置後，再回到原始 TDOA 方程式進行修正。

第 $`i`$ 組 TDOA 量測的殘差：

```math
e_i(p)
=
\lVert p-s_i\rVert_2
-
\lVert p-s_{\mathrm{ref}}\rVert_2
-
\rho_i
```

展開後：

```math
\begin{aligned}
e_i(p)
={}&
\sqrt{
(x-x_i)^2+
(y-y_i)^2+
(z-z_i)^2
}\\
&-
\sqrt{
(x-x_{\mathrm{ref}})^2+
(y-y_{\mathrm{ref}})^2+
(z-z_{\mathrm{ref}})^2
}
-\rho_i
\end{aligned}
```

所有 TDOA 殘差組成：

```math
e(p)=
\begin{bmatrix}
e_1(p)\\
e_2(p)\\
\vdots\\
e_{N-1}(p)
\end{bmatrix}
```

非線性最佳化目標：

```math
\min_p
\frac{1}{2}
\lVert e(p)\rVert_2^2
```

---

## 7. Jacobian Matrix

殘差函數對 $`x`$ 的偏微分：

```math
\frac{\partial e_i}{\partial x}
=
\frac{x-x_i}{\lVert p-s_i\rVert_2}
-
\frac{x-x_{\mathrm{ref}}}
{\lVert p-s_{\mathrm{ref}}\rVert_2}
```

對 $`y`$ 的偏微分：

```math
\frac{\partial e_i}{\partial y}
=
\frac{y-y_i}{\lVert p-s_i\rVert_2}
-
\frac{y-y_{\mathrm{ref}}}
{\lVert p-s_{\mathrm{ref}}\rVert_2}
```

對 $`z`$ 的偏微分：

```math
\frac{\partial e_i}{\partial z}
=
\frac{z-z_i}{\lVert p-s_i\rVert_2}
-
\frac{z-z_{\mathrm{ref}}}
{\lVert p-s_{\mathrm{ref}}\rVert_2}
```

第 $`i`$ 列 Jacobian：

```math
J_i(p)=
\begin{bmatrix}
\dfrac{x-x_i}{\lVert p-s_i\rVert_2}
-
\dfrac{x-x_{\mathrm{ref}}}
{\lVert p-s_{\mathrm{ref}}\rVert_2}
&
\dfrac{y-y_i}{\lVert p-s_i\rVert_2}
-
\dfrac{y-y_{\mathrm{ref}}}
{\lVert p-s_{\mathrm{ref}}\rVert_2}
&
\dfrac{z-z_i}{\lVert p-s_i\rVert_2}
-
\dfrac{z-z_{\mathrm{ref}}}
{\lVert p-s_{\mathrm{ref}}\rVert_2}
\end{bmatrix}
```

---

## 8. Levenberg–Marquardt Method

第 $`k`$ 次迭代的位置為 $`p_k`$，更新量為 $`\Delta p_k`$。

殘差函數的一階近似：

```math
e(p_k+\Delta p_k)
\approx
e(p_k)+J(p_k)\Delta p_k
```

令

```math
e_k=e(p_k)
```

```math
J_k=J(p_k)
```

LM 更新方程式：

```math
(J_k^TJ_k+\lambda I)\Delta p_k
=
-J_k^Te_k
```

求得更新量：

```math
\Delta p_k
=
-(J_k^TJ_k+\lambda I)^{-1}J_k^Te_k
```

實作時直接求解線性系統。

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

更新位置：

```math
p_{k+1}=p_k+\Delta p_k
```

當位置更新量足夠小、殘差不再下降，或達到最大迭代次數時停止。
