#!/bin/bash

if [ ! -r garr.graphml ] ; then
  wget -Ogarr.graphml http://www.topology-zoo.org/files/Garr201201.graphml
fi

