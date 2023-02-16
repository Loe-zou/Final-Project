
#!/bin/bash
for i in *.txt
do
awk '{if (FNR==1) $1=FILENAME; print}' "$i" | sed 's/.txt//g' > "$i.tmp" && mv -f "$i.tmp" "$i"
done

##Credit: eunksung
