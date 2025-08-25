#!/bin/bash

if [ ! -f "./bin/Decompose" ]
then
	make
fi

helpFunction()
{
   echo ""
   echo "Usage: $0 -i INSTANCE -o HIGHEST_ORDER -s SAVE_PATH [-p PRECISION -d]"
   echo -e "\t-i Path of the instance to be decomposed"
   echo -e "\t-o Highest order component to be considered in the decomposition"
   echo -e "\t-s Path in which to save the generated sub-instance folder"
   echo -e "\t-p Precision for the MPFR library (256 bits by default)"
   echo -e "\t-d If set, standardizes the generated sub-instances"
   exit 1 # Exit script after printing help
}

while getopts "i:o:s:p:d" opt
do
   case "$opt" in
      i ) INSTANCE="$OPTARG" ;;
      o ) HIGHEST_ORDER="$OPTARG" ;;
      s ) SAVE_PATH="$OPTARG" ;;
      p ) PRECISION="$OPTARG" ;;
      d ) STANDARD=1 ;;
      ? ) helpFunction ;;
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$INSTANCE" ] || [ -z "$HIGHEST_ORDER" ] || [ -z "$SAVE_PATH" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Default value of optional parameters

if [ -z "$PRECISION" ]
then
   PRECISION=256
fi

if [ -z "$STANDARD" ]
then
   STANDARD=0
fi

shift $((OPTIND -1))

filename=$(basename $INSTANCE)
filename="${filename%.*}"

mkdir -p $SAVE_PATH/$filename

#Decompose input instance
./bin/Decompose $INSTANCE $SAVE_PATH/$filename $HIGHEST_ORDER $PRECISION $STANDARD 0
