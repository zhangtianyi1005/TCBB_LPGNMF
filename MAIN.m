clear
clc
tic
%% Input
load('data/lnc_pro_matrix.mat');%lncRNA-protein interaction matrix
load('data/lnc_matrix.mat');%lncRNA similarity matrix
load('data/pro_matrix.mat');%protein similarity matrix
%%
X=lnc_pro_matrix;
[u1,v1]=size(pro_matrix);
[u2,v2]=size(lnc_matrix);
W1=pro_matrix;
W2=abs(lnc_matrix);

%% obtain the sparse similarity matrices 
 p1=5;
 p2=20;
 sparse_W1=zeros(size(W1));
 sparse_W2=zeros(size(W2)); 
 
for i=1:u1
    row_pro=W1(i,:);
    [tmp,index]=sort(row_pro,'descend');
    p1_nearst=zeros(size(row_pro));
    p1_nearst(index(1:p1))=tmp(1:p1);
    sparse_W1(i,:)=p1_nearst;
end

[u,v]=find(sparse_W1);
for i=1:length(u)
    if(sparse_W1(v(i),u(i))==sparse_W1(u(i),v(i)))
     sparse_W1(v(i),u(i))= sparse_W1(u(i),v(i));
    else
        sparse_W1(u(i),v(i))=0.5*sparse_W1(u(i),v(i));
        sparse_W1(v(i),u(i))= sparse_W1(u(i),v(i));
    end    
end

for i=1:u2
    row_lnc=W2(i,:);
    [tmp,index]=sort(row_lnc,'descend');
    p2_nearst=zeros(size(row_lnc));
    p2_nearst(index(1:p2))=tmp(1:p2);
  
    sparse_W2(i,:)=p2_nearst;
end

[u,v]=find(sparse_W2);
for i=1:length(u)
    if(sparse_W2(v(i),u(i))==sparse_W2(u(i),v(i)))
     sparse_W2(v(i),u(i))= sparse_W2(u(i),v(i));
    else
        sparse_W2(u(i),v(i))=0.5*sparse_W2(u(i),v(i));
        sparse_W2(v(i),u(i))= sparse_W2(u(i),v(i));
    end    
end
%% LPGNMF
rank=10;
lamb1=1;
lamb2=1;
b1=0.25;
b2=0.25;
toc;
[T,W,H]=lpgnmf(X,rank, sparse_W1, sparse_W2,lamb1,lamb2,b1,b2);
result=T;
toc;

%% LOOCV
[pre_label_score]=LOOCV_LPGNMF(X,rank,sparse_W1, sparse_W2,lamb1,lamb2,b1,b2);
%% plot roc
score_predict_bar=pre_label_score(:);
T_ture=X(:);
color0=[0.2,0.5,0.7];
[auc,sn,sp]=roc_curve(score_predict_bar,T_ture,color0);
toc

