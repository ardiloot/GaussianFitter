from __future__ import division, absolute_import, print_function
from numpy.distutils.core import Extension

if __name__ == "__main__":
    from numpy.distutils.core import setup

    ext = Extension(name = "GaussianFitter._Fitter",
                    sources = ["GaussianFitter/src/lmdif.f",
                               "GaussianFitter/src/splev.f",
                               "GaussianFitter/src/gaussian.f90"],
                    extra_compile_args = ["-Wall -O3"],
                    extra_f77_compile_args = ["-Wall -O3"],
                    extra_f90_compile_args = ["-Wall -O3"])
    
    setup(name = "GaussianFitter",
          version = "0.1",
          description = "Uses MINPACK lmdif function for fast fitting of Gaussian functions.",
          author = "Ardi Loot",
          author_email = "ardi.loot@outlook.com",
          url = "",
          packages = ["GaussianFitter",
                      "GaussianFitter.Tests",
                      "GaussianFitter.Examples"],
          ext_modules = [ext],
     )

