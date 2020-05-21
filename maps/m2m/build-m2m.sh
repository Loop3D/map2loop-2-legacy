#!/bin/bash
echo Building map2model .........................................
mkdir /m2m_source/build 
cd /m2m_source/build 
cmake ..
make  
cp map2model ../../../../map2loop/m2m_cpp

