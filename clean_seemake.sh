#!/bin/bash
make clean
find . -iwholename '*cmake*' -not -name CMakeLists.txt -delete
