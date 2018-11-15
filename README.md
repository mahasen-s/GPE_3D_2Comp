# GPE_3D_2Comp
Master branch: 2 Component GPE solver in 3D. Uses adaptively stepped (RK45, Dormand-Prince) symmetrised split-step fourier method. Includes quantum corrections from Bogoliubov theory (Lee-Huang-Yang term).

1comp Branch: 1 Component GPE solver in 3D. Implements higher order solvers (e.g. RK89) and performs these calculations in the interaction picture, which improves convergence. Once the solvers are mature, they will be merged into the multicomponent master branch.
