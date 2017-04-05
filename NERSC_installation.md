# Installing NaMaster at NERSC

The following steps should be followed to install NaMaster at NERSC:

## 1 Install libsharp
libsharp is a C library for spherical harmonic transforms. Follow these steps to install it before moving on to NaMaster:
1. Download libsharp from [its github repository](https://github.com/dagss/libsharp) and unzip the file.
2. From the libsharp folder run `autoreconf-i`, which will generate a `configure` file.
3. Run `./configure --enable-pic` and `make`
4. Create three directories: `bin`, `lib` and `include` in your home directory (unless they're already there). E.g. `mkdir $HOME/bin` etc.
5. Move the contents of `auto/bin`, `auto/lib` and `auto/include` to the corresponding folders you just created in your home directory.

## 2 Edit your .bashrc.ext
From your home directory, open the file `.bashrc.ext` with your favourite text editor and add the following lines at the end of it:
```
export PATH=$HOME/bin:$PATH
export LD_LIBRARY_PATH=$HOME:/usr/common/software/gsl/2.1/intel/lib:/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/lib:$LD_LIBRARY_PATH
export LDFLAGS="-L/usr/common/software/gsl/2.1/intel/lib -L/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/lib -L$HOME/lib"
export CPPFLAGS="-I/usr/common/software/gsl/2.1/intel/include -I/usr/common/software/cfitsio/3.370-reentrant/hsw/intel/include -I$HOME/include"
export CC=cc
```
Note that, if you have previously editted either `$LDFLAGS` or `$CPPFLAGS`, the changes above should be appended to their already existing values to avoid messing up your environment.

## 3 Install NaMaster
Once libsharp has been installed, download or clone NaMaster from its [github repository](https://github.com/damonge/NaMaster) and follow these steps:
1. Run `./configure --prefix=$HOME` and `make`. The last step of this process will fail, but fear not!
2. Run the following command:
```
 cc -std=gnu99 -g -O2 -fopenmp -o .libs/namaster src/namaster-nmt_main.o -L./.libs $LDFLAGS -lnmt -lgsl -lgslcblas -lchealpix -lsharp -lfftpack -lc_utils -lcfitsio -lm
```
3. Run the following commands to manually move things to their corresponding directories:
```
cp .libs/libnmt.* $HOME/lib
cp .libs/namaster $HOME/bin
cp src/namaster.h $HOME/include
```
4. Open the file `setup.py` with your favourite text editor and make sure that the variable `use_icc` is set to `True`.
5. Run `python setup.py install --user`. This will install the python module, `pymaster`.
6. To check that the installation worked, go to the `test` directory and run `python check.py`. If you see a bunch of small numbers and plots coming up after a while (and no errors occurred), you can congratulate yourself: you have a working version of NaMaster on NERSC!
