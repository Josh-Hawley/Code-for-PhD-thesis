%cellPostition.m Is a function that returns a matrix where each row
%corresponds to a single cell, columns correspond to time and the values in
%the matrix gives the cell position the cell is occupying
function posMat=cellPosition(CT,rows,cols)
cells=rows*cols;
for c=1:cells
    [posVec,~]=find(CT==c);
    posMat(c,:)=posVec';
end
posMat=floor((posMat-1)./rows)+1;
end