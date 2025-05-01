j is# Non-Gaussian Surface Simulations

**Author:** Sebastian Korsak

# 🌄 Gaussian Self-Affine Surface Generator

This Python module generates **2D Gaussian random surfaces** with spatial correlations, based on a prescribed **roughness exponent** and **correlation lengths**.  
It implements a spectral method using **FFT-based synthesis** and supports anisotropic scaling.

Adapted from the method described in:

> Yang et al., CMES, vol.103, no.4, pp.251–279, 2014

---

## 🔧 Function

```python
SAimage_fft_2(N, rms, skewness, kurtosis, corlength_x, corlength_y, alpha, non_Gauss=False)
```

---

## 📥 Parameters

| Parameter      | Type    | Description |
|----------------|---------|-------------|
| `N`            | `int`   | Size of the generated surface (`NxN`) |
| `rms`          | `float` | Target **root-mean-square roughness** |
| `corlength_x`  | `float` | Correlation length along the **x-axis** |
| `corlength_y`  | `float` | Correlation length along the **y-axis** |
| `alpha`        | `float` | **Roughness exponent** (also known as Hurst exponent) |
| `non_Gauss`    | `bool`  | Set to `False` to generate **Gaussian surface** only |
| `skewness`, `kurtosis` | `float` | Ignored if `non_Gauss=False` |

---

## 🧠 Mathematical Model
The spatial autocorrelation function is modeled as:

`R(tx, ty) = rms² · exp( - ( sqrt((tx/ξx)² + (ty/ξy)²) )^(2α) )`

Where:
- `ξx`, `ξy`: correlation lengths in x and y directions
- `α`: roughness exponent (e.g., `α = 0.5` for Brownian motion)


---

## 📈 Output

Returns a 2D NumPy array:

```python
surface = SAimage_fft_2(..., non_Gauss=False)
```

- Gaussian-distributed height values
- RMS roughness equal to `rms`
- Zero mean
- Spatial correlations matching `ksix`, `ksiy`, and `alpha`

---

## 🌍 Visualization Example

```python
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

z = SAimage_fft_2(N=256, rms=1.0, skewness=0, kurtosis=3,
                  corlength_x=20, corlength_y=20, alpha=0.8, non_Gauss=False)

X, Y = np.meshgrid(np.arange(z.shape[0]), np.arange(z.shape[1]))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, z, cmap='viridis', edgecolor='none')
plt.title('Gaussian Self-Affine Surface')
plt.show()
```

---

## ⚠️ Notes

- This is a **purely Gaussian** surface generator (zero skewness, kurtosis = 3).
- If you need **non-Gaussian statistics**, set `non_Gauss=True` and provide `skewness` and `kurtosis`. See the full version with rank-order mapping for that.
- Spatial correlation is implemented via FFT filtering based on the **Wiener-Khinchin theorem**.


---

## 📜 References

- Yang, F., et al. *Statistical generation of 3D rough surfaces with arbitrary correlation*, CMES, vol. 103, no. 4, 2014.
- Persson, B.N.J. *Theory of rubber friction and contact mechanics*.

## 👥 Credits

- Python translation by Max Pierini @ [EpiData.it](https://epidata.it)
- Ported from original MATLAB toolbox by Dave (2021)
- Extended and maintained by [your name here]

# Results

![image](https://github.com/user-attachments/assets/be9a6f7b-90c5-4240-89bc-78a295d5c9bb)

![image](https://github.com/user-attachments/assets/a328b67c-366d-4d67-ab96-d95e5e2e6f12)
