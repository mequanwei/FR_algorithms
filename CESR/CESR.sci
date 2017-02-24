///////////////////////////////////////////////////////////////////////////////
// Author:   Ran HE
// Date:     Nov. 2010
// Description:  This code is developed based on Uriel Roque's active set algorithm
//               for the linear least squares problem with nonnegative variables in [1].
//http://www.mathworks.com/matlabcentral/fileexchange/10908-active-set-algorithm/content/activeset.m
//
// Reference:
// [1] Portugal, L.; Judice, J.; and Vicente, L. 1994. A comparison of block
//     pivoting and interior-point algorithms for linear least squares problems
//     with nonnegative variables. Mathematics of Computation 63(208):625-643.
//
// [2] Ran He, Wei-Shi Zheng, Bao-Gang Hu. Maximum Correntropy Criterion for 
//     Robust Face Recognition. IEEE Transactions on Pattern Analysis and Machine
//     Intelligence (TPAMI)£¬in press, 2011.
//
// Copyright (C) 2009-2010 OpenPR
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
//        A - [MxN] matrix 
//        b - [Mx1] vector
//
//  Output:
//        x - solution
//        y - complementary solution
///////////////////////////////////////////////////////////////////////////////

function [x,y] = CESR(A,b)

	[m,n] = size(A);
    weight = ones(m,1);
	F = [];
	G = 1:n;
	x = zeros(n,1);
	Atb = A'*b;
	y = -Atb;

	noready = 1;
	count =1;
	while noready
    
    	//step 1    
    	yG = y(G);
    	if isempty(yG)
        	break;
    	end
    	r = G(yG == min(yG));
     	if(isempty(r))
        	break;
    	end
    	r = r(1);
    
    	if (y(r) < 0) 
        	H1 = [];
        	H2 = r;
        	F = union(setdiff(F,H1),H2);
        	G = union(setdiff(G,H2),H1);
    	else
        	noready = 0;
        	break;
    	end
      

    	noready2 = 1;
    	while noready2;
    
        //step 2
        	AF = A(:,F);
         
         //repmat(sqrt(weight),1, size(AF,2)).*AF;  
         [dt,dt1] = size(AF);
         w1 = sqrt(weight);
         wt1 = AF;
         for i1=1:dt1
           wt1(:,i1)= wt1(:,i1).*w1;
         end
        
         //xF=(wt1)\(weight.*b);
         // disp('here');
         if( rank(wt1'*wt1)<size(wt1'*wt1,1))
           xF = pinv(wt1'*wt1)*(AF'*(weight.*b));
           disp('warning by heran');
         else
           xF = (wt1'*wt1)\(AF'*(weight.*b));
         end
          
         nF = length(xF);
        
        	indt=find(xF < 0);

        	if length(indt)==0  //all(xF >= 0)
        
            	x = [xF; zeros(n-nF,1)];
            	//goto 3
            	noready2 = 0;
            	break;
       
        	else
          
            	index = find(xF<0);
            	t = -x(index)./(xF(index) - x(index));
            	tetha = min(t);
            	r = F(index(t == tetha));
            	if(isempty(r))
            	    break;
            	end
            	r = r(1);
            	x = [(1-tetha)*x(F) + tetha*xF ; zeros(n-nF,1)];
            	H1 = r;
            	H2 = [];
            	F = union(setdiff(F,H1),H2);
            	G = union(setdiff(G,H2),H1);
            	//goto 2
            
        	end
    
  	end //while noready2
   
    
    	//step 3
    	AG = A(:,G);
    	
    [dt,dt1] = size(AG);
        wt1=AG;
    for i1= 1:dt1
      wt1(:,i1)= wt1(:,i1).*weight;   
    end
   
    yG = wt1' * (AF * xF - b);
    y(G) = yG;
    
    //calculate the weight
    weight = (AF *xF-b).^2;
    weight = exp(-weight/(1*mean(weight)));
    
    
    	printf('%d ",count);
    	count=count+1;
    	if(count>200)
        	break;
  	  end
   
	end

	//put results into their corresponding positions
	a = x;
	b = y;
	x = zeros(n,1);
	y = zeros(n,1);
	x(F) = a(1:length(F));
	y(G) = b(1:length(G));
	
endfunction

