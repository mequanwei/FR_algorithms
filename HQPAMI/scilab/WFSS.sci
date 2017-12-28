///////////////////////////////////////////////////////////////////////////////
// Author:   Ran HE (rhe@nlpr.ia.ac.cn)
// Date:     June 2013
//  Description: error detection. A robust sparse representation method based
// on the multiplicative form of half-quadratic minimization. In each iteration, 
// outliers are detected and given small weights. The parameters of 
// different robust estimators will affect robustness such that they should 
// be well tuned for a specific  problem.
// Welsch estimator was currently used. If you want to select another one,
// you can comment the weighting codes around line 136.
//
// This code is developed based on Lee's Efficient sparse coding algorithms in [1].
// http://ai.stanford.edu/~hllee/softwares/nips06-sparsecoding.htm
//
// Reference:
// [1] H. Lee, A. Battle, R. Raina, A. Ng.  Efficient sparse coding algorithms. In
//     Neural information processing systems, 2006,19:801¨C808. 
//
// [2] Ran He, Wei-Shi Zheng, Bao-Gang Hu, XiangWei Kong.
//     A Regularized Correntropy Framework for Robust Pattern Recognition. 
//     Neural Computation (NECO), 2011, 23(8):2074-2100.
//
// [3] Ran He, Wei-Shi Zheng,Tieniu Tan, and Zhenan Sun. 
//     Half-quadratic based Iterative Minimization for Robust Sparse Representation. 
//     IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013. (Accepted)
//
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  Input:
//        A - [MxN] matrix 
//        y - [Mx1] vector
//        gamma1 - the regularization parameter to control sparsity
//        sigmal - the kernel parameter in correntropy to control robustness
//  Output:
//        x - sparse solution
//        weight - auxiliary variable
///////////////////////////////////////////////////////////////////////////////
function [x,weight] = WFSS(A,y,gamma1,sigma1)
AtA = A'*A;
Aty = A'*y;
[L,M] = size(A);

rankA = min(size(A,1)-10, size(A,2)-10);

// Step 1: Initialize
usexinit = 0;

    x= zeros(M,1);
    theta= zeros(M,1);
    act= zeros(M,1);
    allowZero = 0;

fobj = 0; //fobj_featuresign(x, A, y, AtA, Aty, gamma);

ITERMAX=400;
optimality1=0;
weight = ones(L,1);

// A1 = repmat(sqrt(weight),1,M).*A;    
A1=A;
for ti = 1:M
  A1(:,ti) =sqrt(weight).*A1(:,ti),
end

    y1 = sqrt(weight).*y;
//     AtA = A1'*A1;
    Aty = A1'*y1; 

for iter=1:ITERMAX
  // check optimality0
   
    
    act_indx0 = find(act == 0);  
    t_ind = find(abs(x)>0);
    if(length(t_ind)>0)
        grad = A1'*(A1(:,t_ind)*x(t_ind)) - Aty;
    else
        grad =  - Aty;
    end
    theta = sign(x);

    optimality0= 0;
    // Step 2
    [mx,indx] = max (abs(grad(act_indx0)));

    if ~isempty(mx) & (mx >= gamma1) & (iter>1 | ~usexinit) 
        act(act_indx0(indx)) = 1;
        theta(act_indx0(indx)) = -sign(grad(act_indx0(indx)));
        usexinit= 0;
    else
        optimality0= 1;
        if optimality1
            break;
        end
    end
    act_indx1 = find(act == 1);

    if length(act_indx1)>rankA
        warning('sparsity penalty is too small: too many coefficients are activated');
       return;
    end

    if isempty(act_indx1) //length(act_indx1)==0
        if allowZero, allowZero= 0; continue, end
        return;
    end

    
    k=0;

    while 1
        k=k+1;

        if k>ITERMAX
            warning('Maximum number of iteration reached. The solution may not be optimal');
            
            return;
        end

        if isempty(act_indx1) // length(act_indx1)==0
            if allowZero, allowZero= 0; break, end
            return;
        end
        

        // Step 3: feature-sign step
        [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A1, y1, Aty, theta, act, act_indx1, gamma1);
    
//        welsch        
        weight = (A(:,act_indx1)*x(act_indx1)-y).^2;
        delta =sigma1*mean(weight);
        weight = exp(-weight/(delta));
        

                //Huber
//         weight = abs(A(:,act_indx1)*x(act_indx1)-y);
//         thres = 1.5*median(weight);
//         weight(find(weight<=thres))=1;
//         weight(find(weight>thres))=thres./weight(find(weight>thres));
        
        //A1 = repmat(sqrt(weight),1,M).*A;    
        A1=A;
        for ti = 1:M
            A1(:,ti) =sqrt(weight).*A1(:,ti),
        end
        y1 = sqrt(weight).*y;
        Aty = A1'*y1; 
    


        // Step 4: check optimality condition 1

        if optimality1 break; end;
        if lsearch >0 continue; end;

    end

end

if iter >= ITERMAX
    warning('Maximum number of iteration reached. The solution may not be optimal');
end

if 0  // check if optimality
    act_indx1 = find(act==1);
    grad = A1'*(A1*sparse(x)) - Aty;
    norm(grad(act_indx1) + gamma1.*sign(x(act_indx1)))
    find(abs(grad(setdiff(1:M, act_indx1)))>gamma1)
end

fobj = fobj_featuresign(x, A, y, AtA, Aty, gamma1);
endfunction


