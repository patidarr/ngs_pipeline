#!/usr/bin/awk -f

NR<2{print $0;next}{print $0| "sort -r"}
