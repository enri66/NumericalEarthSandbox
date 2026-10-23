#!/bin/bash
cpu=$(top -b -n2 -d2 | grep Cpu | tail -1 | awk '{print $2}')
rss=$(ps -o rss= -C julia | tail -1)
echo "$(hostname -s) us=${cpu} rss_GB=$((rss/1000000))"
