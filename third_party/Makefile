 # This makefile builds the third-party dependencies
 
scip:
	cd res; tar  --keep-old-files -xf scip*; cd scip*; make GMP=true SHARED=true; cd scip*/lib; ln -s libsoplex*.so libsoplex.so; ln -s libzimpl*.so libzimpl.so
	ln -s res/scip*/scip*/ scip

libsbml:
	cd res; unzip libsbml*; cd libsbml*; ./configure --enable-cpp-namespace=libsbml --prefix=`pwd`/build; make; make install
	ln -s res/libsbml*/ libsbml

# mkdir build && cd build
# cmake ..
# make
# make libobjscip
# make libnlpi
# make liblpi
#