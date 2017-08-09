# CNMF_E
Constrained Nonnegative Matrix Factorization for microEndoscopic data 



Download
=======
OPTION 1: download the package using this [LINK](https://github.com/zhoupc/CNMF_E/archive/master.zip)

OPTION 2: (recommended) clone the git repository <https://github.com/zhoupc/CNMF_E.git>. In this way, you are able to get the latest updates of the package within 1-line command or 1-button click. 

Installation
=======
Run cnmfe_setup.m to add CNMF-E package to the search path of MATLAB

`>> cnmfe_setup`

CNMF-E requires [CVX](http://cvxr.com/cvx/) to denoise the extracted calcium traces. CNMF-E will automatically downloaded it and add it to the searching path. 

Other required MATLAB toolboxes are listed below. 

1. /images/images
2. /shared/optimlib/
3. /signal/signal/
4. /stats/stats/
5. /curvefit/curvefit


Example
=======
The best way to get started is running a demo script for analyzing an example data. 

`>> run demos/demo_endoscope.m ` 

Questions
=======
You can ask questions by sending emails to zhoupc1988@gmail.com or joining our [slack channel](https://beat-ica.slack.com) for discussions. 

Reference
=======
**Please cite this paper when you use CNMF-E in your research. Thanks!**

**Zhou, P.**, Resendez, S.L., Rodriguez-Romaguera, J., Jimenez, J.C, Neufeld, S.Q., Stuber, G.D., Hen, R., Kheirbek, M.A., Sabatini, B.L., Kass, R.E., Paninski, L. (2016). [Efficient and accurate extraction of in vivo calcium signals from microendoscopic video data](https://arxiv.org/abs/1605.07266). arXiv Prepr, arXiv1605.07266.

License
=======

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
