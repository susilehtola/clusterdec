#!/bin/bash

# Uncompress recipes
for f in recipes_{2,3,4,5,6,7,8}.dat; do
    if [[ ! -f $f ]]; then
	tar Jxf recipes.tar.xz
    fi
done
