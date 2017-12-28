function [x, theta, act, act_indx1, optimality1, lsearch, fobj] = compute_FS_step (x, A, y, Aty, theta, act, act_indx1, gamma1)

x2 = x(act_indx1);
// A2 = A(:, act_indx1);
AtA2 = A(:,act_indx1)'*A(:,act_indx1);
theta2 = theta(act_indx1);


if(rank(AtA2)<length(AtA2)) then
    yt =( Aty(act_indx1) - gamma1.*theta2 );
    x_new = (AtA2+0.0005*eye(size(AtA2,2)))\yt; 
else
    yt=Aty(act_indx1) - gamma1.*theta2;
    x_new = AtA2\yt;
end

// opts.POSDEF=1; opts.SYM=1; % RR
// x_new = linsolve(AtA2, ( Aty(act_indx1) - gamma1.*theta2 ), opts); % RR
optimality1= 0;
if (sign(x_new) == sign(x2)) then
    optimality1= 1;
    x(act_indx1) = x_new;
    fobj = 0; //fobj_featuresign(x, A, y, AtA, Aty, gamma1);
    lsearch = 1;
    return; 
end

// do line search: x -> x_new
progress = (0 - x2)./(x_new - x2);
lsearch=0;
//a= 0.5*sum((A2*(x_new- x2)).^2);
a= 0.5*sum((A(:, act_indx1)*(x_new- x2)).^2);
b= (x2'*AtA2*(x_new- x2) - (x_new- x2)'*Aty(act_indx1));
fobj_lsearch = gamma1*sum(abs(x2));
[sort_lsearch, ix_lsearch] = sort([progress',1]);
remove_idx=[];
for i = 1:length(sort_lsearch)
    t = sort_lsearch(i); if t<=0 | t>1 continue; end
    s_temp= x2+ (x_new- x2).*t;
    fobj_temp = a*t^2 + b*t + gamma1*sum(abs(s_temp));
    if fobj_temp < fobj_lsearch
        fobj_lsearch = fobj_temp;
        lsearch = t;
        if t<1  remove_idx = [remove_idx ix_lsearch(i)]; end // remove_idx can be more than two..
    elseif fobj_temp > fobj_lsearch
        break;
    else
        if (sum(x2==0)) == 0
            lsearch = t;
            fobj_lsearch = fobj_temp;
            if t<1  remove_idx = [remove_idx ix_lsearch(i)]; end // remove_idx can be more than two..
        end
    end
end


if lsearch >0
    // update x
    x_new = x2 + (x_new - x2).*lsearch;
    x(act_indx1) = x_new;
    theta(act_indx1) = sign(x_new);  // this is not clear...
end

// if x encounters zero along the line search, then remove it from
// active set
if lsearch<1 & lsearch>0
    //remove_idx = find(x(act_indx1)==0);
    remove_idx = find(abs(x(act_indx1)) < 2.2204e-016);
    x(act_indx1(remove_idx))=0;

    theta(act_indx1(remove_idx))=0;
    act(act_indx1(remove_idx))=0;
    act_indx1(remove_idx)=[];
end
fobj_new = 0; //%fobj_featuresign(x, A, y, AtA, Aty, gamma);

fobj = fobj_new;
endfunction
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g] = fobj_featuresign(x, A, y, AtA, Aty, gamma)

f= 0.5*norm(y-A*x)^2;
f= f+ gamma*norm(x,1);

//if nargout >1
    g= AtA*x - Aty;
    g= g+ gamma*sign(x);
//end
endfunction
//%%%%%%%%%%%%%%%%%%%%%

function retval = assert(expr)
retval = true;
    if ~expr 
       // % error('Assertion failed');
        warning ('Assertion failed');
        retval = 0;
    end
endfunction