#!/bin/bash

#Script to identify java VM. Also modifies classpath on cygwin.
#Also used to instantiate headless version of Mesquite.
#Written by Travis Wheeler (www.traviswheeler.com)

#If you find bugs or modify the script to improve functionality, 
# please let us know so we can improve the script for everyone

# Thanks to William Cook (http://www.cs.utexas.edu/users/wcook/) for the bit that
# identifies the path containing this script, and uses it to set the classpath. 
# (the exact site of the work that served as inspiration for this code is lost to antiquity)

#To increase memory allocation adjust the following line:  e.g. to 
#mem="500M" 
mem="2G"

#figure out where java lives 
if [ $NINJA_JAVA_HOME ]
then
  java="$NINJA_JAVA_HOME/bin/java"
elif [ $JAVA_HOME ]
then
  java="$JAVA_HOME/bin/java"
else
  tmp=`java -version 2>&1`
  if echo "$tmp" | grep -q "command not found"  # no "java", so try "jre"
  then
    tmp=`jre -version 2>&1`
    if echo "$tmp" | grep -q "command not found"
    then
       echo "Can't find java. Try setting either the JAVA_HOME environment variable"
       exit
    else
       java="jre"
    fi
  else
   java="java" 
  fi
fi


# figure out where I live, then run java w/ my containing dir as classpath  
dir=`dirname "$0"`
os=`uname`
if test ${os:0:6} = "CYGWIN"
then
  if test ${dir:1:8} = "cygdrive"
  then
    dir="${dir:10:1}:/${dir:12}"
  fi
fi


abspath="$(cd "${0%/*}" 2>/dev/null; echo "$PWD"/"${0##*/}")"
path_only=`dirname "$abspath"`

cmd="$java -server -Xmx$mem -jar $path_only/Ninja.jar $*"
echo "running the following java command:" 1>&2 
echo "$cmd" 1>&2 
$cmd

#When a script ends with an exit that has no parameter, the exit status of the script is the exit status of the last command executed in the script (previous to the exit).
