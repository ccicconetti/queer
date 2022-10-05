#!/bin/bash

if [ "$ROOTDIR" == "" ] ; then
  ROOTDIR=release
fi

base=$(basename $PWD)
main_name=main-$(echo "$base" | cut -f 1 -d _)

if [ -r "$main_name" ] ; then
  read -s -p "The executable file '$main_name' already exists, press any key to overwrite or Ctrl+C to abort" -n 1
  echo
  rm -f "$main_name" 2> /dev/null
fi

target="../../$ROOTDIR/Experiments/$base/$main_name"
if [ ! -x "../../$ROOTDIR/Experiments/$base/$main_name" ] ; then
  echo "could not find the target executable: $target"
else
  ln -s $target $main_name
fi

