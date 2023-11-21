@echo off
rem Copyright 2007-2021 Intel Corporation.
rem 
rem This software and the related documents are Intel copyrighted materials,
rem and your use of them is governed by the express license under which they
rem were provided to you (License). Unless the License provides otherwise,
rem you may not use, modify, copy, publish, distribute, disclose or transmit
rem this software or the related documents without Intel's prior written
rem permission.
rem 
rem This software and the related documents are provided as is, with no
rem express or implied warranties, other than those that are expressly stated
rem in the License.
Rem Intel(R) MPI Library Build Environment

SET I_MPI_ROOT=%~dp0..

title "Intel(R) MPI Library Build Environment for Intel(R) 64 applications"

set library_kind=release

:PARSE_ARGS
if not "%1"=="" (
    if "%1"=="-i_mpi_ofi_internal" (
        set I_MPI_OFI_LIBRARY_INTERNAL=%2
        shift
    )
    if "%1"=="--i_mpi_ofi_internal" (
        set I_MPI_OFI_LIBRARY_INTERNAL=%2
        shift
    )
    if "%1"=="-i_mpi_library_kind" (
        if /i "%2"=="release" (
            set library_kind=%2
        )
        if /i "%2"=="debug" (
            set library_kind=%2
        )
        shift
    )
    if "%1"=="--i_mpi_library_kind" (
        if /i "%2"=="release" (
            set library_kind=%2
        )
        if /i "%2"=="debug" (
            set library_kind=%2
        )
        shift
    )
    if "%1"=="-h" (
        goto :HELP
        shift
    )
    if "%1"=="--help" (
        goto :HELP
        shift
    )
    shift
    goto :PARSE_ARGS
)
goto :EXPORTS

:HELP
echo.
echo "Usage: vars.bat [-i_mpi_ofi_internal [0|1]] [-i_mpi_library_kind [debug|release]]"
echo.
echo "-i_mpi_ofi_internal specifies whether to use libfabric from the Intel(R) MPI Library."
echo.
goto :EOF

:EXPORTS
SET PATH=%I_MPI_ROOT%\bin\%library_kind%;%I_MPI_ROOT%\bin;%PATH%
SET LIB=%I_MPI_ROOT%\lib\%library_kind%;%I_MPI_ROOT%\lib;%LIB%
SET INCLUDE=%I_MPI_ROOT%\include;%INCLUDE%

if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="0" goto :EOF
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="no" goto :EOF
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="off" goto :EOF
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="disable" goto :EOF
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="" goto SET_LIBFABRIC_PATH
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="1" goto SET_LIBFABRIC_PATH
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="yes" goto SET_LIBFABRIC_PATH
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="on" goto SET_LIBFABRIC_PATH
if /i "%I_MPI_OFI_LIBRARY_INTERNAL%"=="enable" goto SET_LIBFABRIC_PATH

:SET_LIBFABRIC_PATH
set PATH=%I_MPI_ROOT%\libfabric\bin\utils;%I_MPI_ROOT%\libfabric\bin;%PATH%
