#!/bin/bash
scp -r `pwd` middle:~/Projects/
ssh middle './startRun.sh'
ssh middle
