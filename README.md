# Solving the Seven-Equation Model (SEM) with Python

This code solves the seven-equation model (SEM) of fluid dynamics.

## Documentation

The most complete documentation is given by the output of `doxygen`; to view this
documentation, navigate to the `doc/doxygen` directory and run `doxygen`:

```
doxygen Doxyfile
```

The root HTML output file is `doc/doxygen/html/index.html`.

## Code Formatting

This code is formatted using `yapf` (Yet Another Python Formatter). The style
file used is `.style.yapf`. To format a single file:

```
yapf -i /path/to/file.py
```

The `-i` argument specifies the formatting to be applied to the file instead of
just printing out the result to the screen. To apply formatting to the entire
code base, just run the provided formatting script:

```
./format.sh
```

## Linter

`pylint3` (`3` corresponding to Python 3) is used as a linter. The configuration
file is `.pylintrc`. To lint a single file and save in a file `pylint_report`,

```
pylint3 /path/to/file.py > pylint_report
```

To apply it to the entire code base,

```
pylint3 sem_python > pylint_report
```
