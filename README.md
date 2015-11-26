# NaMaster

NaMaster implements the MASTER algorithm to compute
angular power spectra in a cut sky.
See http://arxiv.org/abs/astro-ph/0105302 for more
details.

Dependencies:
 - HEALPix: NaMaster uses this pixelization scheme
            only.
 - GSL: used for linear algebra operations.
 
Usage: NaMaster -<opt-name> <option>
Options:
  -map      -> path to file containing HEALPix map
  -map_2    -> path to file containing 2nd map (optional).
               By default it will use the first map.
  -mask     -> path to file containing mask
  -mask_2   -> path to file containing mask for 2nd map (optional).
               By default it will use the first mask.
  -pol      -> spin-0 (0) or spin-2 (1) input map
               If 1, input file must have 2 maps (Q and U)
  -pol_2    -> spin-0 (0) or spin-2 (1) 2nd input map (optional).
               If 1, second input file must have 2 maps (Q and U)
               By default it will use the value of pol.
  -cl_noise -> path to file containing noise Cl(s). The first column
               should be the multipole order l, and the rest should
               correspond to as many power spectra as there are pairs
               of input maps. E.g. if pol=1 and pol_2=1, the columns
               should be, in this order:
                 l  clEE  clEB  clBE  clBB
               These Cl are the pseudo-Cl of the noise for each
               component.
  -coupling -> path to file containing coupling matrix (optional)
               If non-existing, it will be computed and
               written there. Set to "none" if you don't want to save
               the coupling matrix (default).
  -out      -> output prefix
  -nlb      -> number of ells per bin
  -h        -> this help
