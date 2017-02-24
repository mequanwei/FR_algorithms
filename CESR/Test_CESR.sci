///////////////////////////////////////////////////////////////////////////////
// Author:   Ran HE
// Date:     May. 2011
// Description:  Calculate a robust nonnegative sparse representation by CESR
// Reference:
//   Ran He, Wei-Shi Zheng, Bao-Gang Hu. Maximum Correntropy Criterion for 
//   Robust Face Recognition. IEEE Transactions on Pattern Analysis and Machine
//   Intelligence (TPAMI)£¬in press, 2011.
//
// Copyright (C) 2011-2012 OpenPR
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of OpenPR nor the names of its 
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL HOLDER AND CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  Input:
//       data.bin  - number of samples (1 double)
//                 - image height      (1 double)        
//                 - image width       (1 double)     
//                 - probe image       (height*width double)
//                 - gallery images    (119*8*height*width double)
//
//  Example (Scilab Command):
//  clear
//  clc
//  exec('D:\paper\TSR1\CESR.sci');       //change the path before testing
//  exec('D:\paper\TSR1\Test_CESR.sci');  //change the path before testing
//  Test_CESR()
///////////////////////////////////////////////////////////////////////////////

function Test_CESR
  stacksize(20000000);
  //load data AR database from file
  //change the path before testing
  fd1=mopen('D:\paper\TSR1\data.bin','rb'); //the file 'data.bin' exists in the same folder of the source code
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
  [x,y] = CESR(Gallery,probe60);
  
 // plot the indices correponding to nonzero coefficients 
  disp(find(x>0));
 // plot the nonzero coefficients 
  mprintf('%0.2f ',x(find(x>0)));
  //disp(x(find(x>0))');
  
 // plot the maximum coefficient 
  [vx,ix]=max(x);
    printf('\n maximum cofficient: %f; index: %d; Subject ID: %d \n ',vx,ix,ceil(ix/8));
//plot the sparse solution
  ax=1:length(x);
  plot2d(ax, x);
endfunction
