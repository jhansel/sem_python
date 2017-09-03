# Solving the Seven-Equation Model (SEM) with Python

This code solves the seven-equation model (SEM) of fluid dynamics.

## Documentation

The most complete documentation is given by the output of `doxygen`; to view
this documentation, navigate to the `doc/doxygen` directory and run `doxygen`:

```
doxygen Doxyfile
```

The root HTML output file is `doc/doxygen/html/index.html`.

## Installation

First navigate to the `sem_python` directory. To install normally (you will
probably need root privileges):

```
$ python setup.py install
```

This will install the library in the default location. For instructions on
how to customize the install procedure read the output of:

```
$ python setup.py --help install
```

In addition, there are some other commands:

```
$ python setup.py clean # -> will clean all trash (*.pyc and stuff)
$ python setup.py test  # -> will run the complete test suite
$ python setup.py bench # -> will run the complete benchmark suite
$ python setup.py audit # -> will run pyflakes checker on source code
```

To get a full list of available commands, read the output of:

```
$ python setup.py --help-commands
```

To install in developer mode:

```
$ pip install -e .
```
