import math
import unittest
import numpy as np
import GaussianFitter
from scipy import optimize, interpolate

def Gaussian2DNumpy((xs, ys), amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi):
    x, y = xs - x0, ys - y0
    xP = cosPhi * x - sinPhi * y
    yP = sinPhi * x + cosPhi * y
    res = amp * np.exp(-((xP) ** 2.0) / (2.0 * sigmaX ** 2.0) - \
                            ((yP) ** 2.0) / (2.0 * sigmaY ** 2.0)) + offset
    return res

def SplevFast(x, tck):
    t, c, k = tck
    x = x.ravel()
    y, _ = interpolate._fitpack._spl_(x, 0, t, c, k, 0)
    return y

class TestFitting(unittest.TestCase):

    def SciPyFitFuncSigma(self, (xs, ys), amp, x0, y0, sigmaX, sigmaY, offset):
        cosPhi, sinPhi = math.cos(self.phi), math.sin(self.phi)
        return Gaussian2DNumpy((xs, ys), amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi)
    
    def SciPyFitFuncZ(self, (xs, ys), amp, x0, y0, z0, offset):
        cosPhi, sinPhi = math.cos(self.phi), math.sin(self.phi)
        sigmaX = SplevFast(z0, self.tckSigmaX)
        sigmaY = SplevFast(z0, self.tckSigmaY)
        return Gaussian2DNumpy((xs, ys), amp, x0, y0, sigmaX, sigmaY, offset, cosPhi, sinPhi)

    def _CompareWithCurveFitSigma(self, amp, x0, y0, sigmaX, sigmaY, offset):
        np.random.seed(0)
        
        # Mesh
        xM, yM = np.meshgrid(self.xs, self.ys)
        
        # Generate noisy gauss
        noise = 0.05 * np.random.uniform(size = xM.shape)
        gauss = GaussianFitter.Gaussian2D(xM, yM, amp, x0, y0, sigmaX, sigmaY, offset) + noise
        
        # Fit2dGaussSigma
       
        sol1, pcov1, residuals1 = GaussianFitter.Fit2dGaussSigma(xM, yM, gauss, self.initial)
        fitStds1 = np.sqrt(np.diag(pcov1))
    
        # curve_fit
        sol2, pcov2 = optimize.curve_fit(self.SciPyFitFuncSigma, (xM.ravel(), yM.ravel()), gauss.ravel(), self.initial)
        residuals2 = self.SciPyFitFuncSigma((xM, yM), *sol2) - gauss
        fitStds2 = np.sqrt(np.diag(pcov2))
        
        decimal = 7
        np.testing.assert_almost_equal(sol1, sol2, decimal = decimal)
        np.testing.assert_almost_equal(pcov1, pcov2, decimal = decimal)
        np.testing.assert_almost_equal(residuals1, residuals2, decimal = decimal)
        np.testing.assert_almost_equal(fitStds1, fitStds2, decimal = decimal)
        
    def _CompareWithCurveFitZ(self, amp, x0, y0, sigmaX, sigmaY, offset):
        np.random.seed(0)
        
        # Generate test z-dependence
        self.tckSigmaX = interpolate.splrep(self.zExp, self.sigmaXExp)
        self.tckSigmaY = interpolate.splrep(self.zExp, self.sigmaYExp)
        
        # Mesh
        xM, yM = np.meshgrid(self.xs, self.ys)
        
        # Generate noisy gauss
        noise = 0.05 * np.random.uniform(size = xM.shape)
        gauss = GaussianFitter.Gaussian2D(xM, yM, amp, x0, y0, sigmaX, sigmaY, offset) + noise
        
        # Fit2dGaussSigma
        sol1, pcov1, residuals1 = GaussianFitter.Fit2dGaussZ(xM, yM, gauss, self.initial, self.tckSigmaX, self.tckSigmaY)
        fitStds1 = np.sqrt(np.diag(pcov1))
        
        # curve_fit
        sol2, pcov2 = optimize.curve_fit(self.SciPyFitFuncZ, (xM.ravel(), yM.ravel()), gauss.ravel(), self.initial)
        residuals2 = self.SciPyFitFuncZ((xM, yM), *sol2) - gauss
        fitStds2 = np.sqrt(np.diag(pcov2))
        
        decimal = 7
        np.testing.assert_almost_equal(sol1, sol2, decimal = decimal)
        np.testing.assert_almost_equal(pcov1, pcov2, decimal = decimal)
        np.testing.assert_almost_equal(residuals1, residuals2, decimal = decimal)
        np.testing.assert_almost_equal(fitStds1, fitStds2, decimal = decimal)

    def testSigma1(self):
        self.phi = 0.0
        self.xs = np.linspace(-0.2e-6, 0.2e-6, 3)
        self.ys = np.linspace(-0.2e-6, 0.2e-6, 3)
        self.initial = (0.99, 1e-9, 1e-9, 200e-9, 200e-9, 0.1)
        self._CompareWithCurveFitSigma(1.0, 3e-9, 10e-9, 100e-9, 300e-9, 37.0)
    
    def testSigma2(self):
        self.phi = 0.0
        self.xs = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.ys = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.initial = (0.99, 1e-9, 1e-9, 200e-9, 200e-9, 0.1)
        self._CompareWithCurveFitSigma(10.0, 3e-9, 10e-9, 150e-9, 350e-9, 1.0)

    def testSigma3(self):
        self.phi = 0.0
        self.xs = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.ys = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.initial = (0.99, 1e-9, 1e-9, 50e-9, 200e-9, 0.1)
        self._CompareWithCurveFitSigma(0.9, 10e-9, 10e-9, 150e-9, 350e-9, 1.0)
        
    def testZ1(self):
        self.phi = 0.0
        self.xs = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.ys = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.initial = (0.99, 1e-9, 1e-9, 0.0, 0.1)
        
        self.zExp = [-1e-6, 0.0e-6, 0.5e-6, 1e-6]
        self.sigmaXExp = [300e-9, 200e-9, 150e-9, 100e-9]
        self.sigmaYExp = [140e-9, 200e-9, 250e-9, 400e-9]
        
        self._CompareWithCurveFitZ(0.9, 10e-9, 10e-9, 150e-9, 250e-9, 1.0)

    def testZ2(self):
        self.phi = 0.0
        self.xs = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.ys = np.linspace(-0.2e-6, 0.2e-6, 7)
        self.initial = (5.0, 1e-9, 1e-9, 0.0, 0.1)
        
        self.zExp = [-1e-6, 0.0e-6, 0.5e-6, 1e-6]
        self.sigmaXExp = [300e-9, 200e-9, 150e-9, 100e-9]
        self.sigmaYExp = [140e-9, 200e-9, 250e-9, 400e-9]
        
        self._CompareWithCurveFitZ(10.0, 3e-9, 10e-9, 150e-9, 350e-9, 1.0)    
        
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFitting)
    unittest.TextTestRunner().run(suite)