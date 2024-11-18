#cut -d$'\t' -f6,11,4,10 mane_grch38.txt
awk 'BEGIN{FS=OFS="\t"}{print $4"|"$6"|"$11"|"$10, $6, $11, $4, $10}' mane_grch38.txt > mane_grch38_txlist.tsv
