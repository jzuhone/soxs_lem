## SOXS LEM Scripts

This repository collects a set of scripts to use for LEM simulations of various sources. It will be
updated over time. 

* `cgm_example.py`: simulations of the CGM using various galactic halos

# Installing pyXSIM and SOXS

pyXSIM and SOXS can be installed in a few different ways. The simplest way is via the conda package if
you have the [Anaconda Python Distribution](https://store.continuum.io/cshop/anaconda/):

```
[~]$ conda install -c jzuhone pyxsim soxs
```

The second way to install pyXSIM/SOXS is via pip. pip will attempt to download the dependencies and
install them, if they are not already installed in your Python distribution:

```
[~]$ pip install pyxsim soxs
```

**NOTE: For the time being, if you are using Python 3.10 from within Anaconda, you must use pip to install pyXSIM and SOXS.**

Alternatively, to install into your Python distribution from source:
```
[~]$ git clone https://github.com/lynx-x-ray-observatory/soxs
[~]$ cd soxs
[~]$ python -m pip install .

[~]$ git clone https://github.com/jzuhone/pyxsim
[~]$ cd pyxsim
[~]$ python -m pip install .
```