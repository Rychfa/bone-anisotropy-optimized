#!/usr/bin/env bash

FILENAME="fastcode_presentation"
head -n 5 ${FILENAME}.html > tmp.html
echo "<style>" >> tmp.html
cat css/common.css >> tmp.html
echo "</style>" >> tmp.html
tail -n +9 ${FILENAME}.html |\
    grep -v '\$(".remark-slide-number").after' \
    >> tmp.html

if hash decktape 2>/dev/null; then
    decktape tmp.html ${FILENAME}.pdf
elif hash docker 2>/dev/null; then
    docker run --rm -t -v `pwd`:/slides astefanutti/decktape /slides/tmp.html /slides/${FILENAME}.pdf
elif hash singularity 2>/dev/null; then
    singularity run docker://astefanutti/decktape tmp.html ${FILENAME}.pdf

else
    echo >&2 "This script requires 'decktape', 'docker' or 'singularity' to convert HTML to PDF."
    echo >&2 "None of these are installed. Aborting."
    rm tmp.html
    exit 1
fi

rm tmp.html
