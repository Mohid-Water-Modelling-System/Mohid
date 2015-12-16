#!/bin/bash

if [[ `git status --porcelain` ]]; then
    echo "YES"
else
    echo "NO"
fi

