from setuptools import setup, find_packages
import os
# from setuptools import setup, find_packages

#def is_package(path):
#    return (
#        os.path.isdir(path) and
#        os.path.isfile(os.path.join(path, '__init__.py'))
#        )


#def find_packages(path, base=""):
#    """ Find all packages in path """
#    packages = {}
#    for item in os.listdir(path):
#        dir = os.path.join(path, item)
#        if is_package(dir):
#            if base:
#                module_name = "%(base)s.%(item)s" % vars()
#            else:
#                module_name = item
#            packages[module_name] = dir
#            packages.update(find_packages(dir, module_name))
#    return packages

setup(
    name='smolyak',
    version='0.0.1',
    author='EconForge',
    author_email='slyon@stern.nyu.edu',
    packages=find_packages(),
    # packages=['byumcl.interpolate', 'byumcl.misc', 'byumcl.uhlig', 'byumcl.matstat'],
    url='https://github.com/EconForge/Smolyak.jl',
    description='Smolyak Grid interpolation.',
    long_description=open('README.md').read()
)
