# CNMF_E
Constrained Nonnegative Matrix Factorization for microEndoscopic data. 'E' also suggests 'extension'. It is built on top of [CNMF](https://github.com/epnev/ca_source_extraction) with supports to 1 photon data. 

## Download
OPTION 1: download the package using this [LINK](https://github.com/zhoupc/CNMF_E/archive/master.zip)

OPTION 2: (recommended) clone the git repository <https://github.com/zhoupc/CNMF_E.git>. In this way, you are able to get the latest updates of the package within 1-line command or 1-button click. 
```bash 
git clone --recurse-submodules https://github.com/zhoupc/CNMF_E.git 
```
## Installation
Run cnmfe_setup.m to add CNMF-E package to the search path of MATLAB

`>> cnmfe_setup`

<!-- CNMF-E requires [CVX](http://cvxr.com/cvx/) to denoise the extracted calcium traces. CNMF-E will automatically downloaded it and add it to the searching path.  -->

Required MATLAB toolboxes are listed below. 

1. /images/images
2. /shared/optimlib/
3. /signal/signal/
4. /stats/stats/
5. /curvefit/curvefit

## Examples
The best way to get started is running a demo script for analyzing an example data. We included three demo scripts for you to try. You can modify these demos to process your own datasets. 

1. demo_endoscope.m : this demo is good for exploratory analysis of your data. It gives you a good sense of different stages of CNMF-E pipeline and how different parameter selections influence your final results. 

  `>> run demos/demo_endoscope.m ` 

2. demo_large_data_1p.m : this demo is good for processing large-scale dataset with the minimal manual intervention. It can process small data as well and should be used in most automated analysis. 

  `>> run demos/demo_large_data_1p.m`

3. demo_large_data_2p.m : this demo is the same as demo_large_data_1p.m. It is optimized for processing 2p data. 

  `>> run demos/demo_large_data_2p.m`

## Questions
You can ask questions by sending emails to zhoupc1988@gmail.com or joining our [slack channel](https://beat-ica.slack.com) for discussions. An email request for slack channel invitation is required,  or you can ask your collegues who have joined the slack channel already for an invitation. 

Please read more from the [wiki page](https://github.com/zhoupc/CNMF_E/wiki) 

## Python 

This repository is the native implementation of the CNMF-E model & algorithms using MATLAB. There is also a python implementation of the model in [CaImAn](https://github.com/flatironinstitute/CaImAn). However, these two implementations are not exactly the same in some details and we are still working to make them consistent. For those people who love Python but still want to run CNMF-E using the MATLAB implementation, we provided a Python wrapper (thanks to [Tim Machado](https://github.com/tamachado)) for calling MATLAB version CNMF-E in Python.  You can use the ipython notebook file python_wrapper/ analyze_cnmfe_matlab.ipynb to get started. In the meanwhile, [Thomas Akam](https://github.com/ThomasAkam) provided a solution to load CNMF-E results from python, please check [it](https://github.com/ThomasAkam/CNMF-E_mat2npy) out . 

## Reference
**Please cite these papers when you use CNMF-E in your research. Thanks!**

**Zhou, P.**, Resendez, S.L., Rodriguez-Romaguera, J., Jimenez, J.C, Neufeld, S.Q., Giovannucci, A., Friedrich, J., Pnevmatikakis, E.A., Stuber, Garret D ,  Stuber, G.D., Hen, R., Kheirbek, M.A., Sabatini, B.L., Kass, R.E., Paninski, L. (2018). [Efficient and accurate extraction of in vivo calcium signals from microendoscopic video data](https://elifesciences.org/articles/28728). eLife,  pp.e28728.  (When you use the (Corr+PNR) algorithm in the initialization step or use the ring model for estimating background components)

Pnevmatikakis, E.A., Soudry, D., Gao, Y., Machado, T.A., Merel, J., Pfau, D., Reardon, T., Mu, Y., Lacefield, C., Yang, W., Ahrens, M., Bruno, R.,Jessell, T.M., Peterka, D.S., Yuste, R., and Paninski L., 2016. [Simultaneous denoising, deconvolution, and demixing of calcium imaging data](http://www.sciencedirect.com/science/article/pii/S0896627315010843). Neuron, 89(2), pp.285-299. (The original CNMF framework paper; When you use the vanilla CNMF model in this paper ) 



## License
Copyright 2016 Pengcheng Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
