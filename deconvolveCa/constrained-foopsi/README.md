# constrained-foopsi
Implementation of the constrained deconvolution spike inference algorithm in Matlab.

 spike inference using a constrained foopsi approach:   <br />
      min      | sum(sp)    <br />
    c,sp,b,c1  |            <br />
      subject to:  sp|  >= 0     <br />
                  b| >= 0         <br />
                   G*c| = sp        <br />
                    c1| >= 0          <br />
                  ||y-b-c - c_in|| | <= sn*sqrt(T)   <br />


File constained_foopsi.m
   Variables: |   <br />
---------------|-----------------
   y:    |  raw fluorescence data (vector of length(T))     <br />
   c:    |  denoised calcium concentration (Tx1 vector)     <br />
   b:    |  baseline concentration (scalar)                   <br />
  c1:    |  initial concentration (scalar)                    <br />
   g:    |  discrete time constant(s) (scalar or 2x1 vector) <br />
  sn:    |  noise standard deviation (scalar)                 <br />
  sp:    |   spike vector (Tx1 vector)                        <br />

   USAGE:  <br />
   [c,b,c1,g,sn,sp] = constrained_foopsi(y,b,c1,g,sn,OPTIONS)     <br />
   The parameters b,cin,g,sn can be given or else are estimated from the data

   OPTIONS: (stuct for specifying options)   <br />
         p: order for AR model, used when g is not given (default 2)   <br />
    method: methods for performing spike inference  <br />
   available methods: 'dual' uses dual ascent  <br />
                       'cvx' uses the cvx package available from cvxr.com (default)  <br />
                      'lars' uses the least regression algorithm   <br />
                     'spgl1' uses the spgl1 package available from math.ucdavis.edu/~mpf/spgl1/  (usually fastest)  <br />
   bas_nonneg:   flag for setting the baseline lower bound. if 1, then b >= 0 else b >= min(y)   <br />
   noise_range:  frequency range over which the noise power is estimated. Default [Fs/4,Fs/2]  <br />
   noise_method: method to average the PSD in order to obtain a robust noise level estimate  <br />
   lags:         number of extra autocovariance lags to be considered when estimating the time constants  <br />
   resparse: number of times that the solution is resparsened (default 0). Currently available only with methods 'cvx', 'spgl'  <br />
   

The noise is estimated with a power spectral density approach and the time constants from the signal autocovariance. 

The algorithm can also handle missing data (appearing as NaNs in y) due to motion artifacts or for super-resolution approaches

The algorithm is presented in more detail in

Pnevmatikakis, E. A., Gao, Y., Soudry, D., Pfau, D., Lacefield, C., Poskanzer, K., ... & Paninski, L. (2014). A structured matrix factorization framework for large scale calcium imaging data analysis. arXiv preprint arXiv:1409.2903. http://arxiv.org/abs/1409.2903

Toolbox Dependencies
=======
The signal processing toolbox is optional. If present it is used for computing the power spectral density and the autocovariance function.

License
=======

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
