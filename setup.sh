export ISR_UNFOLD_WD=`pwd`

export LD_LIBRARY_PATH=$ISR_UNFOLD_WD/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=$ISR_UNFOLD_WD/TUnfold:$ROOT_INCLUDE_PATH

# TODO use path variable in Linkdef.h

#check if ./lib exist if not create
if [ ! -d "lib" ]; then
  # Control will enter here if $DIRECTORY doesn't exist.
  mkdir lib
fi
