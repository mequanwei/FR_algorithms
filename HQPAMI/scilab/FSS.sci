///////////////////////////////////////////////////////////////////////////////
// Description:  This code is developed by Lee's Efficient sparse coding algorithms in [1].
// http://ai.stanford.edu/~hllee/softwares/nips06-sparsecoding.htm
//
// Reference:
// [1] H. Lee, A. Battle, R. Raina, A. Ng.  Efficient sparse coding algorithms. In
// Neural information processing systems, 2006,19:801¨C808. 

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//  Input:
//        A - [MxN] matrix 
//        y - [Mx1] vector
//        gamma1 - the regularization parameter to control sparsity
//  Output:
//        x - sparse solution
///////////////////////////////////////////////////////////////////////////////

function [x] = FSS(A,y,gamma)
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

ITERMAX=500;
en=[];
optimality1=0;
for iter=1:ITERMAX
    // check optimality0
    act_indx0 = find(act == 0);
    grad = AtA*sparse(x) - Aty;
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
            // if ~assert(max(abs(x))==0), save(fname_debug, 'A', 'y', 'gamma', 'xinit'); error('error'); end
            if allowZero, allowZero= 0; break, end
            return;
        end

        // Step 3: feature-sign step
        [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A, y, Aty, theta, act, act_indx1, gamma);
        //function [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A, y, Aty, theta, act, act_indx1, gamma1)
        
        
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
    grad = AtA*sparse(x) - Aty;
    norm(grad(act_indx1) + gamma.*sign(x(act_indx1)),'inf')
    find(abs(grad(setdiff(1:M, act_indx1)))>gamma)
end

fobj = fobj_featuresign(x, A, y, AtA, Aty, gamma);

