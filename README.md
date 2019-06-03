# Polyquant CT Reconstruction Toolbox
This Matlab toolbox allows direct quantitative reconstruction from polyenergetic X-ray computed tomography (CT) measurements. We hope you find it useful and welcome any feedback or questions [j.mason@ed.ac.uk].

## Features
- Allows quantitative reconstruction into **electron density**, mass density, proton stopping power, quasi-monoenergetic and more.
- Iterative statistical reconstruction under Poisson noise.
- Metal artefact compensation and correction.
- Designed for use with [Michigan Image Reconstruction Toolbox](https://web.eecs.umich.edu/~fessler/code) operators: ensure the toolbox is in your path (by running its 'setup.m') before running these demos.
- Non-negative total variation (TV) reguarisation in 2D and 3D, adapted from [UNLocBoX](https://epfl-lts2.github.io/unlocbox-html/).
- Integrated polyenergetic scatter estimation and mitigation.
- Accelerated ordered sub-sets algorithm with bit-reversal ordering.

## Demos
We have included several demos to cover some of its functionality, including:
- 2D fanbeam CT reconstruction of brain, head, chest, abdomen and pelvis regions.
- 2D fanbeam CT metal artefact mitigation, from double titanium hip implants.
- 3D cone-beam CT reconstruction of head and pelvis regions, under 'full-fan' and 'half-fan' scans respectively.
- 3D cone-beam CT reconstruction with integrated polyenergetic scatter modelling (PolySKS).

## References
These methods are presented in the following publications (please cite if using):
- [Jonathan H Mason et al 2017 Phys. Med. Biol. 62 8739](https://doi.org/10.1088/1361-6560/aa9162)
- [Jonathan H Mason et al 2018 Phys. Med. Biol. 63 225001](https://doi.org/10.1088/1361-6560/aae794)

For more details, extensions and its use in radiotherapy, you can read the thesis:
[Quantitative cone-beam computed tomography reconstruction for radiotherapy planning](http://hdl.handle.net/1842/33193 )

## Acknowledgements
Thanks to Mike Davies, Bill Nailon and Alessandro Perelli for their collaboration and supervision during the development of this work. Another thanks to Alessandro for kindly reviewing this code.
