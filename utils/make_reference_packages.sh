#! /usr/bin/env bash

display_usage() {
  echo "wrong number of arguments, usage:"
  echo "make_reference_packages.sh <references-dir-PhyloMagnet> <output-dir-rpkg>"
}

if [  $# -ne 2 ]
then
		display_usage
		exit 1
fi

mkdir -p $2
cwd=$(pwd)

for d in $1/*/;
do
    id=${d%/}
    cd $d; md5sum * > ${id##*/}.md5; cd $cwd
    tar -czf $2/${id##*/}.tgz -C $d .
done
