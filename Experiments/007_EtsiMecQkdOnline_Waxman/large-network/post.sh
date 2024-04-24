#!/bin/bash

if [ ! -d post ] ; then
  mkdir post
fi

( ../main-007 --print-header ; cat data/*.csv ) | KEYS=algorithm,edge-nodes python3 post.py
