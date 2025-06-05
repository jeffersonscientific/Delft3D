FROM mcr.microsoft.com/windows/servercore:ltsc2022

# Set powershell as shell
SHELL ["powershell", "-Command", "$ErrorActionPreference = 'Stop'; $ProgressPreference = 'Continue'; $verbosePreference='Continue';"]

# Download and install Python 3.11.0
ADD https://www.python.org/ftp/python/3.11.0/python-3.11.0-amd64.exe C:\\python-3.11.0-amd64.exe
RUN Start-Process C:\\python-3.11.0-amd64.exe -ArgumentList '/quiet InstallAllUsers=1 PrependPath=1' -Wait
RUN Remove-Item C:\\python-3.11.0-amd64.exe
RUN New-Item -ItemType SymbolicLink -Path 'C:\\Program Files\\Python311\\python3.exe' -Target 'C:\\Program Files\\Python311\\python.exe'

# Install CMake
ADD https://github.com/Kitware/CMake/releases/download/v3.30.1/cmake-3.30.1-windows-x86_64.msi C:\\cmake-3.30.1-windows-x86_64.msi
RUN Start-Process msiexec -Wait -ArgumentList '/i C:\\cmake-3.30.1-windows-x86_64.msi /qn'
RUN setx /M PATH $($Env:PATH + ';C:\\Program Files\\CMake\\bin')
RUN Remove-Item C:\\cmake-3.30.1-windows-x86_64.msi

# Install Visual C++ Redistributable for Visual Studio 2015-2022 (includes vcruntime140.dll)
ADD https://aka.ms/vs/17/release/vc_redist.x64.exe C:\\vc_redist.x64.exe
RUN Start-Process C:\\vc_redist.x64.exe -ArgumentList '/install /quiet /norestart' -Wait
RUN Remove-Item C:\\vc_redist.x64.exe