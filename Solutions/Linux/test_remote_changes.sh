#!/bin/sh

LOCAL=$(git rev-parse @)
REMOTE=$(git rev-parse @{u})
BASE=$(git merge-base @ @{u})

if [ $LOCAL = $REMOTE ]; then
    echo "Up-to-date"
elif [ $LOCAL = $BASE ]; then
    echo "Ahead"
elif [ $REMOTE = $BASE ]; then
    echo "Behind"
else
    echo "Diverged"
fi
