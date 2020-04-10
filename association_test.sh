plink --bfile $1 --assoc --out assoc_result
sort -k 8 -nr assoc_result.assoc | head -n 20 > top20snps.txt