#!/bin/bash
cp -r ../guarana /tmp/guarana_build  
DIR=/tmp/guarana_build
cd $DIR


# remove code which should not be published
rm -r gap/supple
rm -r tst/completeOldCode
rm -r tst/repsWerner

# remove tilde files 
find . -iname \*~ | xargs -n 20 rm -f

# remove CVS stuff
find . -iname .hg | xargs -n 20 rm -rf

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

tar cf Guarana-$VERS.tar guarana_build
gzip -9 Guarana-$VERS.tar
tar cf Guarana-$VERS.tar guarana_build

# assemble all necessary files
mkdir Guarana
mv Guarana-$VERS.tar Guarana
mv Guarana-$VERS.tar.gz Guarana
cp guarana_build/README Guarana
cp guarana_build/PackageInfo.g Guarana
mv guarana_build Guarana/guarana

#scp -r /tmp/Guarana/* jjm@campbell.mcs.st-andrews.ac.uk:~/public_html/software/Guarana

rsync --delete --verbose --progress --stats --compress --rsh=/usr/bin/ssh --recursive --times --perms --links /tmp/Guarana/* jjm@campbell.mcs.st-andrews.ac.uk:~/public_html/software/Guarana/

# remove the unnecessary diretory
rm -r /tmp/Guarana









