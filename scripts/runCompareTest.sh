#!/bin/bash

# $1 = TransferApp
# $2 = CompareApp
# $3 = Input file
# $4 = Test directory
# $5 = Comparison tolerance
# #6... Additional params for transfer app

TRANSFER=$1
COMPARE=$2
INPUT=$3
TEST_DIR=$4
TOL=$5
shift 5

echo "Transfer app: $TRANSFER"
echo "Compare app: $COMPARE"
echo "Input file: $INPUT"
echo "Test directory: $TEST_DIR"
echo "Tolerance: $TOL"

ps=`xml_grep 'patchfile' $TEST_DIR/$INPUT --text_only`
pfile1=`echo $ps | awk -F ' ' '{print $1}'`
pfile2=`echo $ps | awk -F ' ' '{print $2}'`
rs=`xml_grep 'resultfile' $TEST_DIR/$INPUT --text_only`
rfile1=`echo $rs | awk -F ' ' '{print $1}'`
rfile2=`echo $rs | awk -F ' ' '{print $2}'`
bs=`xml_grep 'boundaryfile' $TEST_DIR/$INPUT --text_only`
bfile1=`echo $bs | awk -F ' ' '{print $1}'`
bfile2=`echo $bs | awk -F ' ' '{print $2}'`

tmpdir=`mktemp -d -t simraXXXXXX`
cp $TEST_DIR/$INPUT $tmpdir
cp $TEST_DIR/$pfile1 $tmpdir
cp $TEST_DIR/$pfile2 $tmpdir
cp $TEST_DIR/$rfile1 $tmpdir
cp $TEST_DIR/$bfile1 $tmpdir

pushd $tmpdir
$TRANSFER $INPUT $@
test $? -eq 0 || exit 1
echo $COMPARE $rfile1 $rfile2 $bfile1 $bfile2 $pfile2 -tol $TOL
$COMPARE $rfile1 $rfile2 $bfile1 $bfile2 $pfile2 -tol $TOL
test $? -eq 0 || exit 1
popd
rm -Rf $tmpdir
