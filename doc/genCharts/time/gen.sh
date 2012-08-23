#!/bin/bash

CURRDIR=$(pwd)

TMPDIR="$CURRDIR/_tmp"
IMGDIR="$CURRDIR/_img"

DATEXT="dat"
IMGEXT="png"

ASSETSDIR="$CURRDIR/assets"
PYDIR="$CURRDIR/../../../"

MADXFILETWISS="$ASSETSDIR/test.twiss.madx"
MADXFILEFORT="$ASSETSDIR/test.fort.madx"
PYFILE="$ASSETSDIR/test.py"

MADX="/afs/cern.ch/group/si/slap/bin/madx"

if [ ! -f $MADX ]; then
  echo "mad-x doesn't exist in your system in the path $MADX"
  exit
fi

if [ -d $TMPDIR ]; then
  rm -r $TMPDIR
fi

if [ -d $IMGDIR ]; then
  rm -r $IMGDIR
fi

mkdir $TMPDIR
mkdir $IMGDIR

cd $TMPDIR
echo "!! This script gives different results for each machine depending on their performance.
However they should be proportional."

for O in $(seq 1 6); do
  echo -n "[gen] Order $O"

  TIMEMAD=$(/usr/bin/time -f "%e"  $MADX < $MADXFILETWISS 2>&1 1> /dev/null)
  echo -n "."
  TIMEPTC=$({ time -p sed -e "s/no=[0-9]/no=$O/" $MADXFILEFORT | $MADX; } 2>&1 | grep real | cut -d " " -f 2)
  echo -n "."

  TIMEFORT=`echo $($PYFILE $PYDIR -f $O | tail -n 1) + $TIMEPTC | bc`
  echo -n "."
  TIMETWISS=`echo $($PYFILE $PYDIR -t $O | tail -n 1) + $TIMEMAD | bc`
  echo -n "."

  echo $O $TIMETWISS $TIMEFORT >> "time.$DATEXT"
  echo ""
done

gnuplot <<EOF
  set term $IMGEXT enhanced
  set output "$IMGDIR/time.$IMGEXT"
  set logscale y
  set ylabel "Time (s)"
  set xlabel "Order"
  p "$TMPDIR/time.$DATEXT" using 1:2 title 'Time twiss' with linespoints,\
    "$TMPDIR/time.$DATEXT" using 1:3 title 'Time fort' with linespoints
EOF
