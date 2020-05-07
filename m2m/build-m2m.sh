cd m2m_source
mkdir build && cd build
cmake ..
make
cp map2model ../../
cd ..
rm -rf build
