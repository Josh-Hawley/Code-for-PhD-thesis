% vecTOmat.m Josh's implementation of the Matlb function vec2mat, because 
% he didn't have the correct toolbox installed at the time! Converts a
% vector to a matrix of user defined column number.

function M=vecTOmat(vector,cols)
rows=length(vector)/cols;
if floor(rows)~=rows
    error('Cols (number of columns) value needs to be a factor of the number of elements in the vector!')
end
M=zeros(rows,cols);
for i=1:cols  
    M(:,i)=vector((i-1)*rows+1:i*rows);
end
end