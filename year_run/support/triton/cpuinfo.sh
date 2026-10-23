#!/bin/bash
echo "$(hostname -s): $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2) | fma flag: $(grep -m1 flags /proc/cpuinfo | tr ' ' '\n' | grep -c -x fma) | avx2: $(grep -m1 flags /proc/cpuinfo | tr ' ' '\n' | grep -c -x avx2) | avx: $(grep -m1 flags /proc/cpuinfo | tr ' ' '\n' | grep -c -x avx)"
