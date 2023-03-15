%randColourMapBright.m Generates a random colormap of size n colours. 
%colourBias should be a vector that biases [r g b] values where 0 is no 
%bias and 1 pushes the values up to a value of 1.

%Example: colourBias=[1 0 0] would make the values appear more red.

function map=randColourMapBright(n,colourBias)

map=rand(n,3);
map(map<0.5)=map(map<0.5)*2;
map(:,3)=map(:,3)*0.8;

map(:,1)=map(:,1) + (1-map(:,1))*colourBias(1);
map(:,2)=map(:,2) + (1-map(:,2))*colourBias(2);
map(:,3)=map(:,3) + (1-map(:,3))*colourBias(3);

end