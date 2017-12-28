%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Ran HE (rhe@nlpr.ia.ac.cn)
% Date:     June. 2013
% Description:  Calculate a robust sparse representation by two forms of
% half-quadratic minimization
% FSS: a fast algorithm to solve Lasso problem
% WFSS: error detection. A robust sparse representation method based on the multiplicative form
% EFSS: error correction. A robust sparse representation method based on the additive form
% Reference:
%[1]Ran He, Wei-Shi Zheng,Tieniu Tan, and Zhenan Sun. 
%   Half-quadratic based Iterative Minimization for Robust Sparse Representation. 
%   IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013. (Accepted)
%[2]Ran He, Wei-Shi Zheng, Bao-Gang Hu, XiangWei Kong.
%   A Regularized Correntropy Framework for Robust Pattern Recognition. 
%   Neural Computation (NECO), 2011, 23(8):2074-2100.
%%

function Test
  load('data');
  
 [x] = FSS(Gallery,probe60,0.005);   
%  % plot the maximum coefficient 
%   [vx,ix]=max(abs(x));
%   sprintf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));

  subplot(1,3,1)
  ax=1:length(x);
  plot(ax, x);
  legend('FSS');

   [x,y] = WFSS(Gallery,probe60,0.005);   
%  % plot the maximum coefficient 
%   [vx,ix]=max(abs(x));
%   sprintf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));

  subplot(1,3,2)
  ax=1:length(x);
  plot(ax, x);
  legend('Multiplicative form');
  
  [x,y] = EFSS(Gallery,probe60,0.005);   
%  % plot the maximum coefficient 
%   [vx,ix]=max(abs(x));
%   sprintf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));

  subplot(1,3,3)
  ax=1:length(x);
  plot(ax, x);
  
  legend('Additive form');



