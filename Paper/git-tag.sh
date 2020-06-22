#!/bin/bash

tag0=__tag.tex
tag1=tag.tex

gitLOG=`git log -1 .`;
gitLOGshort=`git log -1 --abbrev-commit .`;
commitID=`echo "$gitLOG" | grep commit | sed -e"s/commit //g"`;
commitDate=`echo "$gitLOG" | grep Date | sed -e"s/Date:  "//g`;
commitIDshort=`echo "$gitLOGshort" | grep commit | sed -e"s/commit //g"`;

echo "\\def\\commitID{commitID: $commitID}" > $tag0
echo "\\def\\commitDATE{$commitDate}" >> $tag0
echo "\\def\\commitIDshort{commitID: $commitIDshort}" >> $tag0

diffcmd="git diff HEAD -- .";
 diff=`$diffcmd`;
if test -n "$diff"; then
    echo "\\def\\commitSTATUS{UNCLEAN}" >> $tag0
else
    echo "\\def\\commitSTATUS{CLEAN}" >> $tag0
fi


if diff $tag0 $tag1  >/dev/null 2>&1; then
    ## files are identical, no update required
    exit 0;
else
    ## files diff: update
    cp $tag0 $tag1;
fi

