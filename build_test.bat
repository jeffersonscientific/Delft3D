"C:\Program Files\CMake\bin\cmake.exe" ./src/cmake -D CONFIGURATION_TYPE:STRING=TESTS -B build_tests

cd build_tests
"C:\Program Files\CMake\bin\cmake.exe" --build . -j --config Debug
"C:\Program Files\CMake\bin\cmake.exe" --install . --config Debug
cd ..