make install
cp $AMBERHOME/AmberTools/src/pymsmt/tests/gaussian $AMBERHOME/AmberTools/test/pymsmt/mcpb/
cd $AMBERHOME/AmberTools/test/pymsmt/mcpb/
./Run.pymsmt
cp $AMBERHOME/AmberTools/src/pymsmt/tests/gamess $AMBERHOME/AmberTools/test/pymsmt/mcpb/
cd $AMBERHOME/AmberTools/test/pymsmt/mcpb/
./Run.pymsmt
cd $AMBERHOME/AmberTools/src/pymsmt/
