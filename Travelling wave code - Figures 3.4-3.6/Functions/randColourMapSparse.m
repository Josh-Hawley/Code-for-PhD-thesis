%randColourMapSparse.m Produces a colour map where random values 

function map=randColourMapSparse(n,proportion)

rdm=rand(n,1);
rdm(rdm<proportion)=1;
rdm(rdm~=1)=0;

map=[rdm, rdm, rdm];


end