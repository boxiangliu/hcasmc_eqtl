RASQUALDIR=/srv/persistent/bliu2/tools/rasqual/
cd $RASQUALDIR/src
# Not run!  Please export your environment.
export CFLAGS="-I/srv/persistent/bliu2/tools/CLAPACK-3.2.1/INCLUDE -I/srv/persistent/bliu2/tools/CLAPACK-3.2.1/F2CLIBS -I/srv/persistent/bliu2/tools/gsl-2.3/gsl"
export LDFLAGS="-L/srv/persistent/bliu2/tools/CLAPACK-3.2.1 -L/srv/persistent/bliu2/tools/CLAPACK-3.2.1/F2CLIBS -I/srv/persistent/bliu2/tools/gsl-2.3/lib"
make
make install