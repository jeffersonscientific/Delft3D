# GCC Compile Modifications: Raw notes

## On Sherlock:
Rolling code mods below into a new gcc-dev branch, then take it from there...
I have been running setup like, 
$ module use /home/users/myoder96/local/modulefiles/
$ module load spack-sdss-delft3d/
# cmake should now be included in spack-sdss-delft3d/ , but just in case...
$ module load cmake/
$ ./build.sh delft3d4
but that might be because I am using an older tag/branch. The newer version might require a `--compiler` option, then we'll need to be sure we use the correct name -- as we define it in the cmake file. (?)

## Fortran-foo

### Character Arrays, character length issue
Maybe that seg-fault is starting to look like:
/home/groups/bprogers/myoder96/Delft3D/src/tools_gpl/waq_tools/waqpb/waqpb_lib/waqpb_base_settings.f90:46:49:

```
   46 |       procedure(get_accepted_flag_args_interface), deferred :: get_accepted_flag_args
      |                                                 1
Error: Non-polymorphic passed-object dummy argument of ‘get_accepted_flag_args_interface’ at (1)
[
```

Also, see this repo, where SH801 has compiled with GNU in an Ubuntu container and is working on ARM64 compile solution:
https://github.com/SH801/Delft3D-for-Arm64/tree/main

Swan down. Now there's this:

Char Length Error:
```
/home/groups/bprogers/myoder96/Delft3D/src/utils_lgpl/deltares_common/packages/deltares_common/src/unit_utils.f90:77:65:


   77 |           t_unit_category('velocity', ['m/s  ', 'm s-1', 'ms-1 ', 'meter per second']), &
      |                                                                 1
Error: Different CHARACTER lengths (5/16) in array constructor at (1)
```

And it's what it sounds like. I've added the first three Vals out. If we pad them further, to 16 char, this will compile. (one of) the problem(s) with this is that the shorter strings will probably not be recognized correctly at runtime. there is a different way to declare the structure that might eliminate the problem.


### Compile sequence erros:
This (see below) build step (for `swan`, which is included in dfd4-suite, etc.)  appears to run out of sequence, by default. this seems to be pretty common; the build rules are not super smart -- they sometimes just build and link everything in a directory. This also means you can't have a {nm}_bkp.F file, or some other old-school change management -- since it will result in multiple definitions of various things.

gfortran -DCOMMIT_VERSION="" -DHAVE_CONFIG_H=1 -DHAVE_MPI -DSWANEXE -DUSE_MPI -I/home/groups/bprogers/myoder96/Delft3D/src/third_party_open/swan/swan_mpi/../include -I/home/groups/bprogers/myoder96/Delft3D/src/cmake/../version_includes -I/home/groups/bprogers/myoder96/Delft3D/build_swan/fortran_module_dir -I/home/groups/sh_s-dss/share/sdss/spack_sdss/spack/opt/spack/linux-centos7-x86_64/gcc-12.4.0/netcdf-fortran-4.6.1-7beqg4ttloknr6ohvzufoiarq47xvllf/include -I/home/groups/sh_s-dss/share/sdss/spack_sdss/spack/opt/spack/linux-centos7-x86_64/gcc-12.4.0/netcdf-c-4.9.2-bphnug65zxegeb4dxasmrpmjmv6sulna/include -I/home/groups/sh_s-dss/share/sdss/spack_sdss/spack/opt/spack/linux-centos7-x86_64/gcc-12.4.0/mpich-4.2.3-44ww2oq7gnxlxzeb5a7lyd3mk5wprj7m/include -I/home/groups/sh_s-dss/share/sdss/spack_sdss/spack/var/spack/environments/delft3d/.spack-env/view/include -O2 -fPIC -fopenmp -ffixed-line-length-132 -ffree-line-length-512 -fallow-argument-mismatch -cpp -O3 -fPIE -w -c /home/groups/bprogers/myoder96/Delft3D/src/third_party_open/swan/src/serv_xnl4v5.f90 -o CMakeFiles/swan_mpi.dir/home/groups/bprogers/myoder96/Delft3D/src/third_party_open/swan/src/serv_xnl4v5.f90.o

I'll get an `error: file does not exist` (or something....). If I nav to that path and execute the compile command, it will work, then when I re-run `build.sh`, we don't see that error -- because `
serv_xnl4v5.f90.o`  
compiles in the next step. How do we add that as a dependency (in make) or just change it to use the source code (serv_xnl45v.f90) ?

Generally, lots of build sequencing errors -- some of them difficult to trace back. Basically, we need lots of dependency definitions, eg. 
module.o: module.F other_module.o
    Fortran {flags...} module.F other_module.o -o module.o

### Fixed line width compiling:
ok. And here's another thing we're chasing. Looks like this was originally compiled in "Fixed line" mode (`gfortran -ffixed`) and "comments" were just placed off the end of a line, and n_col > col_width, without a comment symbol (!). Of course, not all files are written with the same fixed line length, and the recommended fixed line length arguments, eg from `gnu.cmake`,
`
-ffixed-line-length-132 -ffree-line-length-512
`
are not sufficient (in this example, -ffixed-line-length is much too large). So we see things like this:
        IF ( XP(3,NPO) .GT. 0. ) THEN                                     41.62
           XP(7,NPO) = HS / XP(3,NPO)                                     41.62



that need to be like:
        IF ( XP(3,NPO) .GT. 0. ) THEN                                   !  41.62
           XP(7,NPO) = HS / XP(3,NPO)                                   !  41.62


you can use `vim` to do this, 

:%s/. 41.62/!  41.62/gc
but you have to be careful about replacing the wrong thing, and it takes time...

or you can use this awk-foo to pad all lines and terminate the functional width with a "!"
awk '{printf "%-72s\!%s\n", substr($0,1,72), substr($0, 73)}' /home/groups/bprogers/myoder96/Delft3D/src/third_party_open/swan/src/ocpids.F

BEST way! I think, uses sed:
sed -i 's/^\(.\{72\}\)/\1!/' your_file.txt

test without the `-I` option. Looks like this is maybe the best case. It will add the ! to lines >72 columns, but skip the short lines -- rather than padding them out.

### OMP and MPI:
Basically, some components need to be OMP, some _can_ be OMP, some need/can be MPI, some are probably fully SPP, but there is no one way to compile it, so the CMake modifications need to be at the module level. Here's what I regurgitate from my FD notes:

#### OMP Branch errors:
But then, we get the:
Error: invalid branch to/from OpenMP structured block

all over the place. This is "fixed" by removing the `-fopenmp` flag, so it compiles without OMP capabilities, just MPI. We fix some other errors by forcing the MPI compilers (eg, CC=$(which mpicc) in build.sh ), but this appears to be a legitimate OMP blocking error that will require code modifications to fix.


#### Status and next steps:
LOTS of formatting issues. In particular, comments written off the end of the expected fixed line width -- typically 72-74 characters. Unfortunately, this fixed width is not well defined. 72 is occasionally too short; 74 was too long. Maybe 73? But then other files have lines >>74 and should be compiled with -none. So used `sed` to fix these files (see below), but it's not easy, since you can't just `sed` all the files. You have to identify other ones with(out) the width limit  and treat appropriately.
MPI was not properly referencing in some compile steps. Threw a bunch of CMake-foo at CMakeLists.txt until that worked. I think it's (mostly correct) at this point.
Now getting loads of these -- just in the `swan` module:
/home/groups/bprogers/myoder96/Delft3D/src/third_party_open/swan/src/swancom1.F:2200:132:

 2200 |            IF (STPNOW()) RETURN                                           !40.30
      |                                                                                                                                    ^
Error: invalid branch to/from OpenMP structured block

This look like a legitimate Fortran error. In particular, this looks like an illegal exit from an OMP loop, implemented only if compiled with MPI. This is probably something that the Intel compiler allowed (but probably should not have). It may be a limitation that the code can -- without further reengineering, be built as MPI *or* OMP, not both. 

It appears like this:
!$OMP BARRIER                                                             !40.31
!$OMP MASTER                                                              !40.31
      ARR(1) = HSMN2                                                      !40.30
      ARR(2) = SMN2                                                       !40.30
      CALL SWREDUCE( ARR, 2, SWREAL, SWSUM )                              !40.30
      CALL SWREDUCE( NINDX, 1, SWINT, SWSUM )                             !40.41
#ifdef USE_MPI
      IF (STPNOW()) RETURN                                                !40.30
#endif
      HSMN2 = ARR(1) / REAL(NINDX)                                        !40.30 30.82
       SMN2 = ARR(2) / REAL(NINDX)                                        !40.30 30.82
!$OMP END MASTER                                                          !40.31
!$OMP BARRIER                                                             !40.31

ie, wrapped in an OMP block. My guess is that the whole code block needs to be written separately for OMP and MPI.

I'll report this up to the delft3d developers and see what their plans are.



