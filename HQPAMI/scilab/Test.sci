//////////////////////////////////////////////////////////////////////////////
// Author:   Ran HE (rhe@nlpr.ia.ac.cn)
// Date:     June. 2013
// Description:  Calculate a robust sparse representation by two forms of
// half-quadratic minimization
// FSS: a fast algorithm to solve Lasso problem
// WFSS: error detection. A robust sparse representation method based on the multiplicative form
// EFSS: error correction. A robust sparse representation method based on the additive form
// Reference:
//[1]Ran He, Wei-Shi Zheng,Tieniu Tan, and Zhenan Sun. 
//   Half-quadratic based Iterative Minimization for Robust Sparse Representation. 
//   IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013. (Accepted)
//[2]Ran He, Wei-Shi Zheng, Bao-Gang Hu, XiangWei Kong.
//   A Regularized Correntropy Framework for Robust Pattern Recognition. 
//   Neural Computation (NECO), 2011, 23(8):2074-2100.
////

///////////////////////////////////////////////////////////////////////////////
//  Input:
//       data.bin  - number of samples (1 double)
//                 - image height      (1 double)        
//                 - image width       (1 double)     
//                 - probe image       (height*width double)
//                 - gallery images    (119*8*height*width double)
//
//  Example (Scilab Command):
//clear
//clc
//exec('D:\paper\scilab\HQPAMI\scilab\compute_FS_step.sci');       //change the path before testing
//exec('D:\paper\scilab\HQPAMI\scilab\EFSS.sci');       //change the path before testing
//exec('D:\paper\scilab\HQPAMI\scilab\WFSS.sci');       //change the path before testing
//exec('D:\paper\scilab\HQPAMI\scilab\FSS.sci');
//exec('D:\paper\scilab\HQPAMI\scilab\Test.sci');  //change the path before testing
//Test()
///////////////////////////////////////////////////////////////////////////////

function Test
  stacksize(20000000);
  //load data AR database from file
  //change the path before testing
  fd1=mopen('D:\paper\scilab\HQPAMI\scilab\data.bin','rb'); //the file 'data.bin' exists in the same folder of the source code
  num=mget(1,'d',fd1);    // number of samples
  ir=mget(1,'d',fd1);     // image height
  ic=mget(1,'d',fd1);     // image width
  
  probe60 = mget(ir*ic,'d',fd1);  // read the probe image belonging to the 60th subject
  probe60= probe60';
  
  // read the 119*8 gallery images
  num = num-1;
  Gallery =[];
  for i = 1:num
    s1 = mget(ir*ic,'d',fd1);
    s1= s1';  
    Gallery=[Gallery s1];
  end  
  
  mclose(fd1);
  //disp(size(Gallery));
    
  mprintf('img width: %d; img height: %d \nIteration: ',ic,ir);
 
 [x] = FSS(Gallery,probe60,0.005);   
 // plot the maximum coefficient 
  [vx,ix]=max(abs(x));
  printf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));
//plot the sparse solution
  subplot(1,3,1)
  ax=1:length(x);
  plot2d(ax, x);

   [x,y] = WFSS(Gallery,probe60,0.005,0.5);   
 // plot the maximum coefficient 
  [vx,ix]=max(abs(x));
  printf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));
//plot the sparse solution
  subplot(1,3,2)
  ax=1:length(x);
  plot2d(ax, x);

    [x,y] = EFSS(Gallery,probe60,0.005);   
 // plot the maximum coefficient 
  [vx,ix]=max(abs(x));
  printf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));
//plot the sparse solution
  subplot(1,3,3)
  ax=1:length(x);
  plot2d(ax, x);
endfunction


