#!/bin/bash

# cd pdbs
files=`ls *.PDB`

#dssr=`./dssr`

for file in $files; do
    basename=`basename -s .PDB ${file}`
    ./dssr --input=${file} --json -o=${basename}.json --prefix=temp --more --non-pair --u-turn --po4 --idstr=long 
done

rm -f merged.json
jsons=(`ls *.json`)
last_json=${jsons[${#jsons[@]}-1]}
echo '[' > merged.json
for json in "${jsons[@]}"; do
    cat $json >> merged.json
    if [[ $json != ${last_json} ]]; then
    echo ',' >> merged.json
    fi
done
rm temp-*
echo ']' >> merged.json


exit 0


