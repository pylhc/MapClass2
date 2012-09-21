#!/bin/bash

CURRDIR=$(pwd)

TMPDIR="$CURRDIR/_tmp"
IMGDIR="$CURRDIR/_imgs"

DATEXT="dat"
IMGFORM="png"

ASSETSDIR="$CURRDIR/assets"
PYDIR="$CURRDIR/../../../"

MADXFILE="$ASSETSDIR/test.madx"
PYFILE="$ASSETSDIR/test.py"

MADX="/afs/cern.ch/group/si/slap/bin/madx"

ELEMS=(dr q di s o)
LMIN=0.1
LMAX=7
LINC=0.2

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

echo "Creating graphs for $LMIN < L < $LMAX."

for E in ${ELEMS[@]}; do
  KW=$(egrep "^$E\s*?:" $MADXFILE | cut -d: -f2 | cut -d, -f1 | tr -d '[:space:]')
  echo -n "[gen] Element $KW"
  for L in $(seq $LMIN $LINC $LMAX); do
    echo -n "."
    A=$(echo "$L * 5" | bc)
    sed -e "s/L=/L=$L/" -e "s/ANGLE=15.0/ANGLE=$A/" -e "s/LINE:=/LINE:=\($E\)/" $MADXFILE | $MADX > /dev/null

    CHI2=$($PYFILE $PYDIR | tail -n 1)
    echo $L $CHI2 >> "$E.$DATEXT"
  done
  echo ""

  echo "[gnuplot] $IMGDIR/$E.$IMGFORM"
  gnuplot <<EOF
  # f (x) = a*x**b
  # a = 0; b = 2;
  # fit f(x) '$TMPDIR/$E.$DATEXT' via a,b
  set logscale y
  set format y "10^{%T}"
  set term $IMGFORM enhanced
  set output "$IMGDIR/$E.$IMGFORM"
  set title "$KW ($E)"
  set ylabel "{/Symbol c}^2"
  set xlabel "L"
  p "$TMPDIR/$E.$DATEXT" title 'Data'
  # p "$TMPDIR/$E.$DATEXT" title 'Data', a*x**b title sprintf("Approximation %dx^{%d}",a,b)
EOF
done
