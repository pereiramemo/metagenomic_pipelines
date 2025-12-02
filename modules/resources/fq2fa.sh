awk 'NR % 4 == 1 { sub("@",">",$0); print $0}; NR % 4 == 2 { print $0}' $1
