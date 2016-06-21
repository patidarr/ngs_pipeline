#!/usr/bin/awk -f

NR<1{print $0;next}{print $0| "sort -r"}
