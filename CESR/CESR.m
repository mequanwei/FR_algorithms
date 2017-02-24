
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
    while noready2
        
        
        AF = A(:,F);
        
        [dt,dt1] = size(AF);
        w1 = sqrt(weight);
        wt1 = AF;
        for i1=1:dt1
            wt1(:,i1)= wt1(:,i1).*w1;
        end
        
        if( rank(wt1'*wt1)<size(wt1'*wt1,1))
            xF = pinv(wt1'*wt1)*(AF'*(weight.*b));
            disp('warning by heran');
        else
            xF = (wt1'*wt1)\(AF'*(weight.*b));
        end
        
        nF = length(xF);
        
        indt=find(xF < 0);
        
        if length(indt)==0
            x = [xF; zeros(n-nF,1)];
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
        end    
    end
    AG = A(:,G);
    [dt,dt1] = size(AG);
    wt1=AG;
    for i1= 1:dt1
        wt1(:,i1)= wt1(:,i1).*weight;
    end
    
    yG = wt1' * (AF * xF - b);
    y(G) = yG;  
    weight = (AF *xF-b).^2;
    weight = exp(-weight/(1*mean(weight)));
    count=count+1;
    if(count>200)
        break;
    end
    
end
a = x;
b = y;
x = zeros(n,1);
y = zeros(n,1);
x(F) = a(1:length(F));
y(G) = b(1:length(G));

end

