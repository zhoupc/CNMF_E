# OASIS: Fast online deconvolution of calcium imaging data

Code accompanying the paper "[Fast online deconvolution of calcium imaging data](http://journals.plos.org/ploscompbiol/article?rev=2&id=10.1371/journal.pcbi.1005423)". 

Python: https://github.com/j-friedrich/OASIS

MATLAB: https://github.com/zhoupc/OASIS_matlab

[A brief summary of the FOOPSI approach for denoising & deconvolving calcium imaging](https://github.com/zhoupc/OASIS_matlab/blob/master/document/FOOPSI_.md)
## Get Started
### Installation 
add OASIS function to the search path of MATLAB

`>> oasis_setup`

### Process yoru data 
There is also a high-level function **deconvolveCa.m** for the ease of calling different methods. You only need to specify parameters and then denoise & deconvolve your raw trace. For example 
```matlab 
[c, s, options] = deconvolveCa(y, 'foopsi', 'ar1', 'smin', -3, ...'optimize_pars', true, 'optimize_b', true)
```
In this example, we deconvolve the raw trace $y$ using FOOPSI model and constrain the spike size to be $3\times $ noise levels. The AR coefficients and the baseline were updated automatically. 

For more options, check the **examples/** folder and see the comments in *deconvolveCa.m*. 

### Reproduce figures in the paper
You can reproduce the figures in the paper [3] using the following command (replace * with figure index).   

`>> run examples/paper/fig*.m ` 

## Copyright

Pengcheng Zhou @ Colubmia University 
2018 