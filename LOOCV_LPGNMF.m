function [pre_label_score] = LOOCV_LPGNMF(A,rank,W1,W2,r1,r2,b1,b2)
A_ori = A;
[score_ori] = lpgnmf(A,rank,W1, W2,r1,r2,b1,b2);

index = find(1 == A_ori);
disp(length(index))
for i = 1:length(index)
    A(index(i)) = 0;
    [score] = lpgnmf(A,rank,W1, W2,r1,r2,b1,b2);
    score_ori(index(i)) = score(index(i));
    A = A_ori;
    disp(i)
end
pre_label_score = score_ori;

