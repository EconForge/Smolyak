import os
import numpy

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext



module = 'smoly_cyutils'

setup(cmdclass={'build_ext': build_ext},
      name=module,
      version='1.0',
      ext_modules=[Extension(module, [module + ".pyx"])
                   ],
      include_dirs=[numpy.get_include(),
                    os.path.join(numpy.get_include(), 'numpy')]
      )
