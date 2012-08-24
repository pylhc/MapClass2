#!/bin/bash

CURRDIR=$(pwd)

TMPDIR="$CURRDIR/_tmp"
IMGDIR="$CURRDIR/_imgs"

DATEXT="dat"
IMGFORM="png"

ASSETSDIR="$CURRDIR/assets"
PYDIR="$CURRDIR/../../../"

PYFILE="$ASSETSDIR/errPlot.py"

if [ -d $TMPDIR ]; then
  rm -r $TMPDIR
fi

if [ -d $IMGDIR ]; then
  rm -r $IMGDIR
fi

mkdir $TMPDIR
mkdir $IMGDIR

python $PYFILE $PYDIR

gnuplot <<EOF
  set term $IMGFORM enhanced
  set output "$IMGDIR/dsigma.png"
  set logscale y
  set format y "10^{%T}"
  set ylabel "{/Symbol Ds}/{/Symbol s}_{PTC}"
  set xlabel "Order"
  p "$TMPDIR/sigma.$DATEXT" using 1:2 title "x" with linespoints,\
    "$TMPDIR/sigma.$DATEXT" using 1:3 title "px" with linespoints,\
    "$TMPDIR/sigma.$DATEXT" using 1:4 title "y" with linespoints,\
    "$TMPDIR/sigma.$DATEXT" using 1:5 title "py" with linespoints
EOF
