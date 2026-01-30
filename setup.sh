#!/bin/bash

# Create environment from file (handles everything)

conda env create -f nsforest-nf/container/nsforest/context/nsforest.yml

#or to update an existing environment

conda env update -f nsforest-nf/container/nsforest/context/nsforest.yml
echo ""
echo "Setup complete!"
echo "To activate a new environment, run:"
echo "  conda activate nsforest"
echo " otherwise you are good to go"
