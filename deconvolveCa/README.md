# OASIS: Fast online deconvolution of calcium imaging data

Code accompanying the paper "Fast Active Set Method for Online Spike Inference from Calcium Imaging". [arXiv, 2016]

Python: https://github.com/j-friedrich/OASIS

MATLAB: https://github.com/zhoupc/OASIS_matlab

### Use OASIS 

The implementation of OASIS requires several functions for different formulation of the optimization problems. They can be found within the folder `OASIS_matlab/oasis/`. However, we highly recommend you use the wrapper function `deconvolveCa.m` to run differnet algorithms with one command. 

|                       | AR1  | AR2  | convolution ( difference of 2 expfunctions) | convolution |
| --------------------- | ---- | ---- | ---------------------------------------- | ----------- |
| FOOPSI[1]             | √    | √    | √                                        | √           |
| Constrained-FOOPSI[2] | √    |      |                                          |             |
| Thresholded-FOOPSI[3] | √    | √    |                                          |             |



### Installation
Run setup.m to add OASIS function to the search path of MATLAB

`>> setup`

### Examples
The scripts to reproduce the figures are in the subfolder 'examples' with names obvious from the arXiv paper. They can be run with 

`>> run examples/paper/fig*.m ` 

There is also a function **deconvolveCa.m** wrapping all formulations of the problem. You might only need this one for your problem. See a list of examples  in the demo script **oasis_test.m**, 

`>> edit examples/all_examples.m`

### Benchmarks



### References

[1]Vogelstein, J.T., Packer, A.M., Machado, T.A., Sippy, T., Babadi, B., Yuste, R. and Paninski, L., 2010. Fast nonnegative deconvolution for spike train inference from population calcium imaging. *Journal of neurophysiology*,*104*(6), pp.3691-3704.

[2]Pnevmatikakis, E.A., Soudry, D., Gao, Y., Machado, T.A., Merel, J., Pfau, D., Reardon, T., Mu, Y., Lacefield, C., Yang, W. and Ahrens, M., 2016. Simultaneous denoising, deconvolution, and demixing of calcium imaging data. *Neuron*, *89*(2), pp.285-299.

[3]Friedrich, J., Zhou, P., and Paninski, L., 2016. Fast Active Set Methods for Online Deconvolution of Calcium Imaging Data. *arXiv preprint arXiv:*1609.00639.