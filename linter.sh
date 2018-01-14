#!/bin/sh

outfile=./linter_output

echo "Linting codebase..."

pylint3 sem_python > ${outfile}

echo "Linter output saved to the following file: ${outfile}"
