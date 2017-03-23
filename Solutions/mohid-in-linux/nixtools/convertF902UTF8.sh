#!/bin/bash

# echo colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

to_encode='utf-8'


HELP(){
           echo
           echo " Convert encoding of given files.F90 to utf-8"
           echo
           echo " Usage: $0 [-h] [-d <directory>]"
           echo "    -d|dir                   : Directory to convert"
           echo "    -h|-help|--help          : Show this help"
           echo
}
if [ $# -lt 1 ] ;then
  HELP
  exit 0
fi

if [[ $1 == '-h' || $1 == '-help' || $1 == '--help' ]]; then
    HELP
    exit 0
fi

if [[ $1 == '-d' && ! -z $2 ]]; then
    dir=$2
else
    dir="."
fi

#clear

## get all F90 files in subdirectories
files=$(find $dir -type f -iname "*.F90" -exec file {} \;|grep CRLF|tr ":" " " |cut -d" " -f1)

## 
for file in $files; do
    ## get mime-encoding from each file
    encoding=$(file -b --mime-encoding $file)
    if [[ $encoding == *unknown* ]] ; then
      printf  " convert %-100s    %20s  to  utf-8   ${RED}ERROR${NC} \n" "$file" "$encoding"
      echo
      echo -e " <file -b --mime-encoding> command can't guess the mime-encoding for $file\n Please convert manually to utf-8 and run this script again"
      echo
      exit
    else
      printf  " convert %-100s    %20s  to  utf-8   ${GREEN}OK${NC} \n" "$file" "$encoding"
      iconv  -f $encoding -t utf-8 "$file" -o "$file.utf8"
      dos2unix < "$file" |  iconv  -f $encoding -t utf-8 > "$file.utf8"
    fi
    rm $file
    mv "$file.utf8" $file
done

