import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import johnsonsu
from skimage.filters import threshold_otsu

# Custom function (assumed you have it in the same directory)
from j_johnson_M import f_johnson_M

def SAimage_fft_2(N=500, rms=3, skewness=1, kurtosis=3, corlength_x=10, corlength_y=10, alpha=0.9,
                  threshold=1000.0, non_Gauss=True, corr=False):

    # --- Gaussian Surface Generation --- #
    txmin = tymin = -N/2
    txmax = tymax = N/2
    dtx = dty = (txmax - txmin) / N
    tx = np.arange(txmin, txmax, dtx)
    ty = np.arange(tymin, tymax, dty)

    R = np.zeros((N+1, N+1))
    for i, txx in enumerate(tx):
        for j, tyy in enumerate(ty):
            r = np.sqrt((txx/corlength_x)**2 + (tyy/corlength_y)**2)
            R[i, j] = rms**2 * np.exp(-np.abs(r)**(2*alpha))

    FR = np.fft.fft2(R, s=[N, N])
    AMPR = np.sqrt(dtx**2 + dty**2) * np.abs(FR)

    X = np.random.rand(N, N)
    X = (X - np.mean(X)) / np.std(X)
    XF = np.fft.fft2(X, s=[N, N])

    YF = XF * np.sqrt(AMPR)
    z = np.real(np.fft.ifft2(YF, s=[N, N]))
    z = (z - np.mean(z)) * rms / np.std(z)
    z_gs = z.copy()

    if not non_Gauss:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X, Y = np.meshgrid(np.arange(N), np.arange(N))
        ax.plot_surface(X, Y, z_gs, cmap='viridis', edgecolor='none')
        ax.view_init(30, 150)
        plt.title('Gaussian Surface')
        plt.show()

        plt.imshow(z_gs < threshold, cmap='gray')
        plt.title('Binarized Gaussian Surface')
        plt.show()
        return z_gs - np.min(z_gs)  # Return surface shifted to be positive

    # --- Non-Gaussian Noise --- #
    coef, j_type, err = f_johnson_M(0, rms, skewness, kurtosis)
    gamma, delta, xi, lam = coef
    print(f"JohnsonSU params: gamma={gamma}, delta={delta}, xi={xi}, lambda={lam}")

    z_ngn = johnsonsu.rvs(a=gamma, b=delta, loc=xi, scale=lam, size=(N, N))

    # Clip and normalize the noise
    z_ngn = np.clip(z_ngn, *np.percentile(z_ngn, [0.1, 99.9]))
    z_ngn = (z_ngn - np.mean(z_ngn)) / np.std(z_ngn) * rms

    # --- Rank-based Mapping --- #
    v_gs = z_gs.flatten()
    v_ngn = z_ngn.flatten()

    Igs = np.argsort(v_gs)
    vs_ngn_sorted = np.sort(v_ngn)

    v_ngs = np.zeros_like(v_gs)
    v_ngs[Igs] = vs_ngn_sorted
    z_ngs = v_ngs.reshape((N, N))

    # Mirror and shift to positive values
    z_ngs = -z_ngs
    z_ngs = z_ngs - np.min(z_ngs)

    # --- Plot Final Surface --- #
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(np.arange(N), np.arange(N))
    ax.plot_surface(X, Y, z_ngs, cmap='hot', edgecolor='black')
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)
    ax.set_zlim(np.min(z_ngs), np.max(z_ngs))
    ax.view_init(30, 150)
    plt.title('Non-Gaussian Surface')
    plt.show()

    # --- Binarization --- #
    plt.imshow(z_ngs > threshold, cmap='gray')
    plt.title('Binarized Non-Gaussian Surface')
    plt.show()

    # --- Optional Correlation Function --- #
    if corr:
        hhcf1d = np.zeros(N//2)
        for ndif in range(N//2):
            surf1 = z_ngs[:, :N-ndif]
            surf2 = z_ngs[:, ndif:]
            hhcf1d[ndif] = np.sqrt(np.mean((surf1 - surf2)**2))

        plt.loglog(np.arange(N//2), hhcf1d)
        plt.grid()
        plt.xlabel('log(r(nm))')
        plt.ylabel('log(G(r) (nm))')
        plt.title('1-D height-height correlation function')
        plt.show()
    
    return z_ngs