I. General information

The material model is characterized by the following parameters

  1. Number of a 'periodic piece' - a collection of layers that is repeated over vertical axis.

  2. The range - the maximal interaction radius for potentials used in the experiment

  3. The length of the modeled material piece

  4. The height of the modeled material piece

  5. Atom types

The energy of the modeled piece is defined as a sum of energies of all atoms in this piece. Interactions in the range of the potential are computed.


The piece of material is defined by the following set of parameters:

(h0, d0, w0, ..., h(k-1), d(k-1), w(k-1))

for k layers numbered from 0 to k - 1
hi - the distance between the i-th layer and the i-1'th layer, 
di - displacement of the first atom in the layer, 
wi - interatomic distance in the layer.


II. JSON mappings

{
# Material's model definition
  "model" : {
    #number of layers
    "nlay" : value,

    # Model's range
    "range" : value,

    # Length of the modeled piece
    "length" : value,

    # height of the modeled piece
    "height" : value,
    
    #atom types
    "atoms" : [...]
  }

# Lattice parameters (vector of 3 * number of layers parameters

    "lattice" : [...]

# Search box parameters
    "box" : {
# Smallest point       
       "a" : [...], 
# Biggest point
       "b": [...]
}
}
