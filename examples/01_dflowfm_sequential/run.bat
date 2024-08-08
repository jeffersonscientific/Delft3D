@ echo off

set dimrdir=p:\d-hydro\dimrset\latest

call %dimrdir%\x64\dimr\scripts\run_dimr.bat dimr_config.xml

    rem To prevent the DOS box from disappearing immediately: remove the rem on the following line
rem pause
