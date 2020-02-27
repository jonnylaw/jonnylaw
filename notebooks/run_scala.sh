#!/bin/bash

mkdir -p notebooks/data
for file in notebooks/*.sc
  do
    amm $file
  done
