# `advec_mom` in Halide

This is a rough port of a part of the [CloverLeaf miniapp](http://warwick-pcav.github.io/CloverLeaf/)
to [Halide](http://halide-lang.org).

## Current Caveats
- Need to get the boundaries correct.  Currently, the call to
  CloverLeaf's `advec_mom` and the Halide version don't do the
  exact same computation due to the boundaries.
- Need to ensure correctness
- Need to implement the rest of the cases for this kernel
- Autotuning
