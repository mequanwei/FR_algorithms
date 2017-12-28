//////////////////////////////////////////////////////////////////////////////
// Author:   Ran HE (rhe@nlpr.ia.ac.cn)
// Date:     June 2013
// Description: error correction. A robust sparse representation method based
// on the additive form of half-quadratic minimization. In each iteration, 
// outliers are detected and corrected. The parameters of 
// different robust estimators will affect robustness such that they should 
// be well tuned for a specific problem. Welsch estimator was currently used.
//
// This code is developed based on Lee's Efficient sparse coding algorithms in [1].
// http://ai.stanford.edu/~hllee/softwares/nips06-sparsecoding.htm
//
// Reference:
// [1] H. Lee, A. Battle, R. Raina, A. Ng.  Efficient sparse coding algorithms. In
//     Neural information processing systems, 2006,19:801¨C808. 
//
// [2] Ran He, Wei-Shi Zheng,Tieniu Tan, and Zhenan Sun. 
//     Half-quadratic based Iterative Minimization for Robust Sparse Representation. 
//     IEEE Trans. on Pattern Analysis and Machine Intelligence, 2013. (Accepted)
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//  Input:
//        A - [MxN] matrix 
//        y - [Mx1] vector
//        gamma1 - the regularization parameter to control sparsity
//        sigmal - the kernel parameter in correntropy to control robustness
//  Output:
//        x - sparse solution
//        delta - auxiliary variable for errors
//////////////////////////////////////////////////////////////////////////////

function [x,delta] = EFSS(A,y,gamma)
AtA = A'*A;
Aty = A'*y;
[L,M] = size(A);

l0norm=[];

rankA = min(size(A,1)-10, size(A,2)-10);

// Step 1: Initialize
usexinit = 0;

    x= (zeros(M,1));
    theta= (zeros(M,1));
    act= (zeros(M,1));
    allowZero = 0;




fobj = 0; //fobj_featuresign(x, A, y, AtA, Aty, gamma);

ITERMAX=1200;
optimality1=0;

    A1 = A;    
    y1 = y;
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

    if ~isempty(mx) & (mx >= gamma) & (iter>1 | ~usexinit)
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

    if isempty(act_indx1) %length(act_indx1)==0
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
        [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A1, y1, Aty, theta, act, act_indx1, gamma);
    
        //welsch
        delta = A(:,act_indx1)*x(act_indx1)-y;
        delta1 = delta.^2;
        delta1 = delta1/(0.6*mean(delta1));
        delta = delta.*(1-exp(-delta1));

        //log cosh
        //delta = A(:,act_indx1)*x(act_indx1)-y;
       // delta = delta - 0.5*tanh(0.5*delta);
       
       //L1-L2
//        delta = A(:,act_indx1)*x(act_indx1)-y;
//        delta = delta.*(1-1./sqrt(1+delta.^2/mean(1*delta.^2)));
        
//fair
// delta = A(:,act_indx1)*x(act_indx1)-y;
//         ad = 1*abs(delta)/mean(abs(delta));
//         delta = delta.*(1-1./(0.5+ad));

// Huber
//         delta = A(:,act_indx1)*x(act_indx1)-y;
//         thres= 1.5* median(abs(delta));
//         delta(find(abs(delta)<=thres))=0;
//         delta(find(abs(delta)>thres))= delta(find(abs(delta)>thres)) -sign(delta(find(abs(delta)>thres)))*thres;


//         delta = delta.^2;
//         delta = delta/mean(delta);
//         delta = delta.*(1-1./sqrt(1+delta));
//         delta = delta.*(1-exp(-0.5*delta./mean(delta)));

        y1 = y+delta;
        Aty = A1'*y1; 
    
        l0norm =[l0norm length(find(abs(delta)>0.001))];

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
    norm(grad(act_indx1) + gamma.*sign(x(act_indx1)),'inf')
    find(abs(grad(setdiff(1:M, act_indx1)))>gamma)
end

fobj = fobj_featuresign(x, A, y, AtA, Aty, gamma);

// plot(1:5:length(l0norm(1:5:end))*5,l0norm(1:5:end));


