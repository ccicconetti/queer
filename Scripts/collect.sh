#!/bin/bash

if [ "$SERVERS" == "" ] ; then
  echo "empty environment variable SERVERS"
  exit 1
fi

if [ "$DATADIR" == "" ] ; then
  DATADIR=data
fi

if [ -d $DATADIR ] ; then
  read -n 1 -p "There is already a directory '$DATADIR', do you want to proceed? (Ctrl+C to abort)"
else
  mkdir $DATADIR
fi

if [[ -d $DATADIR && "$CLEANDIR" != "" ]] ; then
  rm -rf $DATADIR/* >& /dev/null
fi

for i in $SERVERS ; do
  scp -r $i:$PWD/$DATADIR $DATADIR-$i
done

for i in $(cd $DATADIR-$(echo $SERVERS | cut -f 1 -d ' ') && ls -1) ; do
  cat $DATADIR-*/$i > $DATADIR/$i
done

for i in $SERVERS ; do
  rm -rf $DATADIR-$i
done
