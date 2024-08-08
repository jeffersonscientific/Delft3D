@echo off
REM List all subfolders in the current directory
for /d %%D in (*) do (
    echo Checking folder: %%D
    cd %%D
    if exist run.bat (
        echo Running run.bat in %%D
        call run.bat
    ) else (
        echo No run script in %%D
    )
    cd ..
)
