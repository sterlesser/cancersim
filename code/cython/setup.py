from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("forcefunccelltypes", ["forcefunccelltypes.pyx"])]

setup(
  name = 'forcefunccelltypes',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)