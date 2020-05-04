#!/bin/sh

cat > GromacsCmdLineArgs.py << EOF
import argparse

def get_commandLineArguments(cmdline):
    """Parse command line arguments of GROMACS
       using argparse

    Returns:
        the dictionary of arguments
    """
    cmdlineparser = argparse.ArgumentParser()
EOF

awk '
NR%2==1{
    option=$1
    if ($2 ~ /\(/ || $2 ~ /\)/){
        default=$2;
        helpline1=$3;
        helpline2=$4 " " $5 " " $6;
    } else if ($2 ~ /\[/ && $3 ~ /\]/) {
        default=$4;
        helpline1=$2 " " $3;
        helpline2=$5 " " $6 " " $7;
    } else {
        if ($3 ~ /\(/ && $3 ~ /\)/){
            helpline1=$2;
            default=$3;
            helpline2=$4 " " $5 " " $6;
        } else if ($3 ~ /\(/ && $4 ~ /\)/){
            helpline1=$2;
            default=$3 " " $4;
            helpline2=$5 " " $6;
        } else if ($3 ~ /\(/ && $5 ~ /\)/){
            helpline1=$2;
            default=$3 " " $4 " " $5;
            helpline2=$6 " " $7;
        } else {
            helpline1=$2;
            default=$3;
            helpline2=$4 " " $5 " " $6;
        }
    }
}
NR%2==0{
        if(option ~ /\[(no)\]/){
	newoption = option;
	sub(/\[(no)\]/, "", newoption);
	printf("%s%s","    cmdlineparser.add_argument('\''",newoption);
        printf("%s%s","'\'', nargs='\''?'\'', const=","True");
	printf("%s%s%s\n",", default='\''",substr(default,2,length(default)-2),"'\'',");
        printf("%s%s%s%s%s\n","        help='\''",helpline1,helpline2,$0,"'\'')");
	newoption = option;
	sub(/\[/, "", newoption);
	sub(/\]/, "", newoption);
	printf("%s%s","    cmdlineparser.add_argument('\''",newoption);
        printf("%s%s","'\'', nargs='\''?'\'', const=","False");
	printf("%s%s%s\n",", default='\''",substr(default,2,length(default)-2),"'\'',");
        printf("%s%s%s%s%s\n","        help='\''",helpline1,helpline2,$0,"'\'')");
    } else {
        printf("%s%s","    cmdlineparser.add_argument('\''",option);
        printf("%s%s","'\'', nargs='\''?'\'', const='\''",substr(default,2,length(default)-2));
        printf("%s%s%s\n","'\'', default='\''",substr(default,2,length(default)-2),"'\'',");
        printf("%s%s%s%s%s\n","        help='\''",helpline1,helpline2,$0,"'\'')");
    }
}' < GromacsCmdLineArgs.txt >> GromacsCmdLineArgs.py

cat >> GromacsCmdLineArgs.py << EOF
    args = cmdlineparser.parse_args(str(cmdline).split())
    return vars(args)

EOF

