import math
import numpy as np
import _Fitter

def Gaussian2D(xs, ys, amp, x0, y0, sigmaX, sigmaY, offset, phi = 0.0):
    # Flatten arrays
    originalShape = None
    if len(xs.shape) != 1:
        originalShape = xs.shape
        xs = xs.ravel()
        ys = ys.ravel()
        
    cosPhi, sinPhi = math.cos(phi), math.sin(phi)
    res = _Fitter.gaussian2d(xs, ys, amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi)
    
    # Restore shape
    if originalShape is not None:
        res = res.reshape(originalShape)
    return res

def Fit2dGaussSigma(xs, ys, values, initial, phi = 0.0, \
                    ftol = 1.49012e-8, xtol = 1.49012e-8, gtol = 0.0, \
                    maxfev = 0, factor = 100):
    # Flatten arrays
    originalShape = None
    if len(xs.shape) != 1:
        originalShape = xs.shape
        xs = xs.ravel()
        ys = ys.ravel()
        values = values.ravel()
        
    # Init parameters
    n = len(initial)
    if maxfev == 0:
        maxfev = 100 * (n + 1)
    
    # Do fitting
    cosPhi, sinPhi = math.cos(phi), math.sin(phi)
    fitRes = _Fitter.fit_gaussian_sigma(xs, ys, values, initial, cosPhi, \
                                        sinPhi, ftol, xtol, gtol, maxfev, factor)
    
    # Process result
    solution, residual, info, nfev, fjac, ipvt, qtf = fitRes  # @UnusedVariable
    shapeFjac = fjac.shape
    fjac = fjac.ravel(order = "F").reshape(shapeFjac[::-1])
    
    # Calculate pcov
    pcov = _CalcPcov(initial, info, fjac, ipvt, residual, values, n)
    
    # Restore shape
    if originalShape is not None:
        residual = residual.reshape(originalShape)
    
    return solution, pcov, residual

def Fit2dGaussZ(xs, ys, values, initial, tckSigmaX, tckSigmaY, phi = 0.0, \
                    ftol = 1.49012e-8, xtol = 1.49012e-8, gtol = 0.0, \
                    maxfev = 0, factor = 100):
    # Flatten arrays
    originalShape = None
    if len(xs.shape) != 1:
        originalShape = xs.shape
        xs = xs.ravel()
        ys = ys.ravel()
        values = values.ravel()
        
    # Init parameters
    n = len(initial)
    if maxfev == 0:
        maxfev = 100 * (n + 1)
    
    # Do fitting
    cosPhi, sinPhi = math.cos(phi), math.sin(phi)
    tSigmaX, cSigmaX, kSigmaX = tckSigmaX
    tSigmaY, cSigmaY, kSigmaY = tckSigmaY
    fitRes = _Fitter.fit_gaussian_z(xs, ys, values, initial, cosPhi, sinPhi, \
                                    tSigmaX, cSigmaX, kSigmaX, 0, \
                                    tSigmaY, cSigmaY, kSigmaY, 0, \
                                    ftol, xtol, gtol, maxfev, factor)
    
    # Process result
    solution, residual, info, nfev, fjac, ipvt, qtf = fitRes  # @UnusedVariable
    shapeFjac = fjac.shape
    fjac = fjac.ravel(order = "F").reshape(shapeFjac[::-1])
    
    # Calculate pcov
    pcov = _CalcPcov(initial, info, fjac, ipvt, residual, values, n)
    
    # Restore shape
    if originalShape is not None:
        residual = residual.reshape(originalShape)
    
    return solution, pcov, residual

def _CalcPcov(initial, info, fjac, ipvt, fvec, values, n):
    pcov = None
    if info in [1, 2, 3, 4]:
        perm = np.take(np.eye(n), ipvt - 1, 0)
        r = np.triu(np.transpose(fjac)[:n, :])
        R = np.dot(r, perm)
        try:
            pcov = np.linalg.inv(np.dot(np.transpose(R), R))
        except (np.linalg.LinAlgError, ValueError):
            pass

        if values.size > len(initial):
            cost = np.sum(fvec ** 2)
            s_sq = cost / (values.size - len(initial))
            pcov = pcov * s_sq
        else:
            pcov.fill(np.inf)
    return pcov
    

if __name__ == "__main__":
    from scipy import interpolate
    import pylab as plt
    
    np.random.seed(0)
    
    # Generate test z-dependence
    zExp = [-1e-6, 0.0e-6, 0.5e-6, 1e-6]
    sigmaXExp = [300e-9, 200e-9, 150e-9, 100e-9]
    sigmaYExp = [140e-9, 200e-9, 250e-9, 400e-9]
    
    tckSigmaX = interpolate.splrep(zExp, sigmaXExp)
    tckSigmaY = interpolate.splrep(zExp, sigmaYExp)

    # Mesh
    xs = np.linspace(-1e-6, 1e-6, 30)
    ys = np.linspace(-1e-6, 1e-6, 31)
    xM, yM = np.meshgrid(xs, ys)
    
    # Generate noisy gauss
    amp, x0, y0, sigmaX, sigmaY, offset = 1.0, 0.0, 0.0, 150e-9, 250e-9, 0.0
    noise = 0.05 * np.random.uniform(size = xM.shape)
    gauss = Gaussian2D(xM, yM, amp, x0, y0, sigmaX, sigmaY, offset) + noise
    
    # Fitting      
    initial = [1.1, 3e-9, 3e-9, 0.0, 0.0]
    #sol1, pcov1, residuals1 = Fit2dGaussSigma(xM, yM, gauss, initial)
    sol1, pcov1, residuals1 = Fit2dGaussZ(xM, yM, gauss, initial, tckSigmaX, tckSigmaY)
    #fitStds1 = np.sqrt(np.diag(pcov1))
    print sol1
    
    # Plotting
    zs = np.linspace(-1e-6, 1e-6, 200)
    sigmaXs = interpolate.splev(zs, tckSigmaX)
    sigmaYs = interpolate.splev(zs, tckSigmaY)

    plt.figure()
    plt.plot(1e6 * zs, 1e9 * sigmaXs)
    plt.plot(1e6 * zs, 1e9 * sigmaYs)
        
    plt.figure()
    plt.pcolormesh(xs, ys, gauss)
        
    plt.show()