#!/bin/bash -eu

docker run --rm -v $(pwd):/data -it hrektts/doxygen /bin/bash -c "doxygen Doxyfile"
