#!/bin/bash
files=`find ~/.mozilla -name "*lock"`
for file in `echo $files`
do
  echo "removing $file..."
  rm "$file"
done
