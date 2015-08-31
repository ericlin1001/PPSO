#!/bin/sh
rm a.txt
make 2>&1 |tee a.txt
#vim a.txt
