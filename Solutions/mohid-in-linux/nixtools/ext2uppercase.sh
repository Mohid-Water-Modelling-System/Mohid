#!/usr/bin/bash

cd ../../../Software

## to lowercase
#ext=F
#find . -type f  -iname "*.${ext}*" -exec sh -c 'a=$(echo "$0" | sed -r "s/([^.]*)\$/\L\1/"); [ "$a" != "$0" ] && mv "$0" "$a" ' {} \;

## to uppercase
ext=f
find . -type f  -iname "*.${ext}*" -exec sh -c 'a=$(echo "$0" | sed -r "s/([^.]*)\$/\U\1/"); [ "$a" != "$0" ] && mv "$0" "$a" ' {} \;

