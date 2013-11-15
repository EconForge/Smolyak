docs will go here

Note that we borrow a theme from readthedocs -- thanks guys.

In order to build the documentation you will need to install the sphinx read the docs plugin. You can do this using

```
pip install sphinx_rtd_theme
```

We have created a new make instruction to automatically take the latest documentation from the master branch, build the html, and publish to github pages so the documentation is available at http://econforge.github.io/Smolyak.jl/index.html. To use this make instruction just enter `make gh-pages` from this directory. When it is finished, you will be returned to the master branch -- all as if nothing happened.
