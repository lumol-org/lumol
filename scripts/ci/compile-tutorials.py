#!/usr/bin/env python
# coding=utf-8
import os
import subprocess
import glob

TUTORIALS = glob.glob('{}/tutorials/*'.format(os.curdir))

if __name__ == '__main__':
    for tutorial in TUTORIALS:
        subprocess.call('Cargo build', shell=True, cwd=tutorial)

