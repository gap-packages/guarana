#!/bin/bash
cp -r ../guarana /tmp  
DIR=/tmp/guarana
cd $DIR


# remove code which should not be published
rm -r gap/supple
rm -r tst/completeOldCode
rm -r tst/repsWerner

# remove tilde files 
find . -iname \*~ | xargs -n 20 rm -f

# remove CVS stuff
find . -iname CVS | xargs -n 20 rm -rf

# remove TODO files
find . -iname TODO | xargs -n 20 rm -f

# remove unnecessary doc files
cd $DIR/doc
rm -f manual.ind manual.dvi manual.log manual.aux manual.idx manual.bbl manual.ilg

#check the version number
cd $DIR
VERS=`cat VERSION`

# remove yourself
cd $DIR
rm pack-guarana.sh 

# create tar archive and compress it
cd /tmp

tar cf Guarana-$VERS.tar guarana
gzip -9 Guarana-$VERS.tar
tar cf Guarana-$VERS.tar guarana

# assemble all necessary files
mkdir Guarana
mv Guarana-$VERS.tar Guarana
mv Guarana-$VERS.tar.gz Guarana
cp guarana/README Guarana
cp guarana/PackageInfo.g Guarana
mv guarana Guarana

#scp -r /tmp/Guarana/* jjm@campbell.mcs.st-andrews.ac.uk:~/public_html/software/Guarana

rsync --delete --verbose --progress --stats --compress --rsh=/usr/bin/ssh --recursive --times --perms --links /tmp/Guarana/* jjm@campbell.mcs.st-andrews.ac.uk:~/public_html/software/Guarana/

# remove the unnecessary diretory
rm -r /tmp/Guarana









