"C:\Program Files\CMake\bin\cmake.exe" ./src/cmake -D CMAKE_BUILD_TYPE=Release -D CONFIGURATION_TYPE:STRING=all -B build_all -D CMAKE_INSTALL_PREFIX=".\build_all\x64\Release\"
 
cd build_all
"C:\Program Files\CMake\bin\cmake.exe" --build . -j --target install --config Release > out.txt
cd ..