export LD_LIBRARY_PATH=/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/lib:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=/home/jhkim/ISR2016/unfolding/TUnfoldISR2016/TUnfold:$ROOT_INCLUDE_PATH

#check if ./lib exist if not create
if [ ! -d "lib" ]; then
  # Control will enter here if $DIRECTORY doesn't exist.
  mkdir lib
fi
