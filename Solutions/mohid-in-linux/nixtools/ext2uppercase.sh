#!/usr/bin/bash

cd ../../../Software
ext=F

## to lowercase
#find . -type f  -iname "*.${ext}*" -exec sh -c 'a=$(echo "$0" | sed -r "s/([^.]*)\$/\L\1/"); [ "$a" != "$0" ] && mv "$0" "$a" ' {} \;

## to uppercase
find . -type f  -iname "*.${ext}*" -exec sh -c 'a=$(echo "$0" | sed -r "s/([^.]*)\$/\U\1/"); [ "$a" != "$0" ] && mv "$0" "$a" ' {} \;

