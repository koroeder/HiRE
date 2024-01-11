#!/bin/bash

DIR="$( dirname -- "${BASH_SOURCE[0]}")"


DIR_ABS="$( realpath -e -- "$DIR")"

echo "Annotating file ${1}"

cut -c 6- ${DIR_ABS}/CurrentOrder.dat > annot.tmp

paste ${1} annot.tmp > ${1}.annotated

rm annot.tmp
