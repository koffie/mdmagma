#!/usr/bin/env bash
ls *.m | parallel --joblog logs/verify_joblog.txt -j ${1:-10} 'magma -n {} >| logs/{.}.txt'
