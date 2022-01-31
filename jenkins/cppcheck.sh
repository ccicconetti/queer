#!/bin/bash

if [ "$TEXT" == "1" ] ; then
  OPTS="--output-file=report.txt"
else
  OPTS="--xml --xml-version=2 --output-file=report.xml"
fi

CPPCHECK=/usr/local/bin/cppcheck
OPTIONS="--library=gnu --library=posix \
         --platform=unix64 --std=c++17 \
         --suppressions-list=cppcheck.suppression \
         --enable=all \
         --suppress=syntaxError \
         -j 8 \
         --force \
         $OPTS"
         
${CPPCHECK} ${OPTIONS} \
  ../QuantumRouting \
  ../Experiments/001_Constant_Rate_PPP/ \
  ../Test \
  ../support/Support \
  ../support/Test
