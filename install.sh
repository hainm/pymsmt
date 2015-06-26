#!/bin/sh

make install
echo "GAUSSIAN TEST..."
cp -r $AMBERHOME/AmberTools/src/pymsmt/tests/gaussian/ $AMBERHOME/AmberTools/test/pymsmt/mcpb/
cd $AMBERHOME/AmberTools/test/pymsmt/mcpb/
./Run.pymsmt
echo "GAMESS-US TEST..."
cp -r $AMBERHOME/AmberTools/src/pymsmt/tests/gamess/ $AMBERHOME/AmberTools/test/pymsmt/mcpb/
cd $AMBERHOME/AmberTools/test/pymsmt/mcpb/
./Run.pymsmt
cd $AMBERHOME/AmberTools/src/pymsmt/
