%%%min||X-W'H||+r1tr(WL1W')+r2tr(HL2H')+b1sum||Wi||+b2sum||Hj||
function [T,W,H]=lpgnmf(X,k,Pw,Ph,r1,r2,b1,b2)
[n,m]=size(X);
a=sum(Pw,2);
b=sum(Ph,2);
Dw=diag(a);
Dh=diag(b);
%% Initialize non-nagetive matrices randomly
W=abs(rand(k,n));  
H=abs(rand(k,m));
%% update process
E=ones(k,k);
W_container=W;
H_container=H;
delta1=1;
delta2=1;


   while  ((delta1>1e-4)|(delta2>1e-4))
      Hup=W*X+r2*H*Ph;
      Hdown=W*W'*H+r2*H*Dh+b2*E*H;
      H=H.*Hup./(Hdown+eps);%%% update H£»
      Wup=H*X'+r1*W*Pw;
      Wdown=H*H'*W+r1*W*Dw+b1*E*W;
      W=W.*Wup./(Wdown+eps);%%% updateW
 
    delta1 = max(sum((W_container - W).^2));
    delta2 = max(sum((H_container - H).^2));
 
    W_container=W;
    H_container=H;
    
   end
   T=W'*H;
end