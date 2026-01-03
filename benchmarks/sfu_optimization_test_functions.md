# Local Optimization Test Functions

A reference for functional forms, gradients, and Hessians of common benchmark functions.

------------------------------------------------------------------------

## 1. Sphere Function

**Functional Form:**

$$f(\mathbf{x}) = \sum_{i=1}^{d} x_i^2$$

**Gradient:**

$$\frac{\partial f}{\partial x_i} = 2x_i$$

$$\nabla f(\mathbf{x}) = 2\mathbf{x}$$

**Hessian:**

$$H = 2I_d$$

where $I_d$ is the $d \times d$ identity matrix.

**Properties:** Global minimum at $\mathbf{x}^* = \mathbf{0}$ with $f^* = 0$. Convex, separable, unimodal.

------------------------------------------------------------------------

## 2. Bohachevsky Function (Function 1)

**Functional Form:** (2 dimensions)

$$f(x_1, x_2) = x_1^2 + 2x_2^2 - 0.3\cos(3\pi x_1) - 0.4\cos(4\pi x_2) + 0.7$$

**Gradient:**

$$\frac{\partial f}{\partial x_1} = 2x_1 + 0.9\pi\sin(3\pi x_1)$$

$$\frac{\partial f}{\partial x_2} = 4x_2 + 1.6\pi\sin(4\pi x_2)$$

**Hessian:**

$$H = \begin{bmatrix} 2 + 2.7\pi^2\cos(3\pi x_1) & 0 \\ 0 & 4 + 6.4\pi^2\cos(4\pi x_2) \end{bmatrix}$$

**Properties:** Global minimum at $\mathbf{x}^* = (0, 0)$ with $f^* = 0$. Multimodal, separable.

------------------------------------------------------------------------

## 3. Sum of Squares Function

**Functional Form:**

$$f(\mathbf{x}) = \sum_{i=1}^{d} i \cdot x_i^2$$

**Gradient:**

$$\frac{\partial f}{\partial x_i} = 2i \cdot x_i$$

$$\nabla f(\mathbf{x}) = \begin{bmatrix} 2x_1 \\ 4x_2 \\ \vdots \\ 2d \cdot x_d \end{bmatrix}$$

**Hessian:**

$$H = \text{diag}(2, 4, 6, \ldots, 2d)$$

$$H_{ij} = \begin{cases} 2i & \text{if } i = j \\ 0 & \text{otherwise} \end{cases}$$

**Properties:** Global minimum at $\mathbf{x}^* = \mathbf{0}$ with $f^* = 0$. Convex, separable, unimodal.

------------------------------------------------------------------------

## 4. Rotated Hyper-Ellipsoid Function

**Functional Form:**

$$f(\mathbf{x}) = \sum_{i=1}^{d}\sum_{j=1}^{i}x_j^2 = \sum_{j=1}^{d}(d - j + 1)x_j^2$$

Expanded: $f(\mathbf{x}) = d \cdot x_1^2 + (d-1) \cdot x_2^2 + \cdots + x_d^2$

**Gradient:**

$$\frac{\partial f}{\partial x_j} = 2(d - j + 1)x_j$$

$$\nabla f(\mathbf{x}) = \begin{bmatrix} 2d \cdot x_1 \\ 2(d-1) \cdot x_2 \\ \vdots \\ 2x_d \end{bmatrix}$$

**Hessian:**

$$H = \text{diag}(2d, 2(d-1), \ldots, 2)$$

$$H_{ij} = \begin{cases} 2(d - i + 1) & \text{if } i = j \\ 0 & \text{otherwise} \end{cases}$$

**Properties:** Global minimum at $\mathbf{x}^* = \mathbf{0}$ with $f^* = 0$. Convex, separable, unimodal, ill-conditioned.

------------------------------------------------------------------------

## 5. Trid Function

**Functional Form:**

$$f(\mathbf{x}) = \sum_{i=1}^{d}(x_i - 1)^2 - \sum_{i=2}^{d}x_i x_{i-1}$$

**Gradient:**

$$\frac{\partial f}{\partial x_1} = 2(x_1 - 1) - x_2$$

$$\frac{\partial f}{\partial x_i} = 2(x_i - 1) - x_{i-1} - x_{i+1} \quad \text{for } 1 < i < d$$

$$\frac{\partial f}{\partial x_d} = 2(x_d - 1) - x_{d-1}$$

**Hessian:**

Tridiagonal matrix:

$$H = \begin{bmatrix} 2 & -1 & 0 & \cdots & 0 \\ -1 & 2 & -1 & \cdots & 0 \\ 0 & -1 & 2 & \cdots & 0 \\ \vdots & & \ddots & \ddots & -1 \\ 0 & \cdots & 0 & -1 & 2 \end{bmatrix}$$

$$H_{ij} = \begin{cases} 2 & \text{if } i = j \\ -1 & \text{if } |i - j| = 1 \\ 0 & \text{otherwise} \end{cases}$$

**Properties:** Global minimum at $x_i^* = i(d + 1 - i)$ with $f^* = -\frac{d(d+4)(d-1)}{6}$. Convex, non-separable, unimodal. Origin is a saddle point.

------------------------------------------------------------------------

## 6. Powell Singular Function

**Functional Form:** (4 dimensions)

$$f(\mathbf{x}) = (x_1 + 10x_2)^2 + 5(x_3 - x_4)^2 + (x_2 - 2x_3)^4 + 10(x_1 - x_4)^4$$

**Gradient:**

$$\frac{\partial f}{\partial x_1} = 2(x_1 + 10x_2) + 40(x_1 - x_4)^3$$

$$\frac{\partial f}{\partial x_2} = 20(x_1 + 10x_2) + 4(x_2 - 2x_3)^3$$

$$\frac{\partial f}{\partial x_3} = 10(x_3 - x_4) - 8(x_2 - 2x_3)^3$$

$$\frac{\partial f}{\partial x_4} = -10(x_3 - x_4) - 40(x_1 - x_4)^3$$

**Hessian:**

Let $u = x_1 - x_4$ and $v = x_2 - 2x_3$.

$$H_{11} = 2 + 120u^2$$

$$H_{12} = H_{21} = 20$$

$$H_{13} = H_{31} = 0$$

$$H_{14} = H_{41} = -120u^2$$

$$H_{22} = 200 + 12v^2$$

$$H_{23} = H_{32} = -24v^2$$

$$H_{24} = H_{42} = 0$$

$$H_{33} = 10 + 48v^2$$

$$H_{34} = H_{43} = -10$$

$$H_{44} = 10 + 120u^2$$

**Properties:** Global minimum at $\mathbf{x}^* = (0, 0, 0, 0)$ with $f^* = 0$. Non-separable. Hessian is singular at the solution.

------------------------------------------------------------------------

## 7. Rosenbrock Function

**Functional Form:**

$$f(\mathbf{x}) = \sum_{i=1}^{d-1}\left[100(x_{i+1} - x_i^2)^2 + (1 - x_i)^2\right]$$

**Gradient:**

$$\frac{\partial f}{\partial x_1} = -400x_1(x_2 - x_1^2) - 2(1 - x_1)$$

$$\frac{\partial f}{\partial x_i} = 200(x_i - x_{i-1}^2) - 400x_i(x_{i+1} - x_i^2) - 2(1 - x_i) \quad \text{for } 1 < i < d$$

$$\frac{\partial f}{\partial x_d} = 200(x_d - x_{d-1}^2)$$

**Hessian:**

The Hessian is a tridiagonal-like sparse matrix. Let $r_i = x_{i+1} - x_i^2$.

Diagonal elements:

$$H_{11} = -400(x_2 - x_1^2) + 800x_1^2 + 2 = -400r_1 + 800x_1^2 + 2$$

$$H_{ii} = 200 + 1200x_i^2 - 400x_{i+1} + 2 = 200 - 400r_i + 800x_i^2 + 2 \quad \text{for } 1 < i < d$$

$$H_{dd} = 200$$

Off-diagonal elements (only adjacent variables are coupled):

$$H_{i,i+1} = H_{i+1,i} = -400x_i \quad \text{for } i = 1, \ldots, d-1$$

All other entries are zero.

**Properties:** Global minimum at $\mathbf{x}^* = (1, 1, \ldots, 1)$ with $f^* = 0$. Non-convex for $d > 2$, non-separable, has a narrow curved valley.

------------------------------------------------------------------------

## 8. Dixon-Price Function

**Functional Form:**

$$f(\mathbf{x}) = (x_1 - 1)^2 + \sum_{i=2}^{d} i(2x_i^2 - x_{i-1})^2$$

**Gradient:**

Let $g_i = 2x_i^2 - x_{i-1}$ for $i \geq 2$.

$$\frac{\partial f}{\partial x_1} = 2(x_1 - 1) - 2 \cdot 2 \cdot g_2 = 2(x_1 - 1) - 4(2x_2^2 - x_1)$$

$$\frac{\partial f}{\partial x_i} = 8i \cdot x_i \cdot g_i - 2(i+1) \cdot g_{i+1} \quad \text{for } 1 < i < d$$

$$= 8i \cdot x_i(2x_i^2 - x_{i-1}) - 2(i+1)(2x_{i+1}^2 - x_i)$$

$$\frac{\partial f}{\partial x_d} = 8d \cdot x_d \cdot g_d = 8d \cdot x_d(2x_d^2 - x_{d-1})$$

**Hessian:**

The Hessian has a banded structure with bandwidth 2.

$$H_{11} = 2 + 4 = 6$$

For $1 < i < d$:

$$H_{ii} = 8i(6x_i^2 - x_{i-1}) + 2(i+1) = 48i \cdot x_i^2 - 8i \cdot x_{i-1} + 2(i+1)$$

$$H_{dd} = 8d(6x_d^2 - x_{d-1}) = 48d \cdot x_d^2 - 8d \cdot x_{d-1}$$

Off-diagonal elements:

$$H_{i-1,i} = H_{i,i-1} = -8i \cdot x_i \quad \text{for } i = 2, \ldots, d$$

$$H_{i,i+1} = H_{i+1,i} = -8(i+1) \cdot x_{i+1} \quad \text{for } i = 1, \ldots, d-1$$

**Properties:** Global minimum at $x_i^* = 2^{-\frac{2^i - 2}{2^i}}$ with $f^* = 0$. Unimodal, non-separable, has a narrow valley.

------------------------------------------------------------------------

## Summary Table

| Function | Dimensions | Minimum Location | Minimum Value | Convex | Separable |
|----|----|----|----|----|----|
| Sphere | $d$ | $\mathbf{0}$ | $0$ | Yes | Yes |
| Bohachevsky | $2$ | $(0, 0)$ | $0$ | No | Yes |
| Sum of Squares | $d$ | $\mathbf{0}$ | $0$ | Yes | Yes |
| Rotated Hyper-Ellipsoid | $d$ | $\mathbf{0}$ | $0$ | Yes | Yes |
| Trid | $d$ | $x_i = i(d+1-i)$ | $-\frac{d(d+4)(d-1)}{6}$ | Yes | No |
| Powell Singular | $4$ | $\mathbf{0}$ | $0$ | No | No |
| Rosenbrock | $d$ | $\mathbf{1}$ | $0$ | No | No |
| Dixon-Price | $d$ | $x_i = 2^{-(2^i-2)/2^i}$ | $0$ | No | No |
