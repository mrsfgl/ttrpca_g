load M:/Documents/MATLAB/mbttda_temp/Data_file/CoilData.mat

var = 0.05;
I = Data + var*randn(size(Data));
sz = size(I);

for i=1:sz(3)
    for j=1:sz(4)
        I(:,:,i,j)=imnoise(I(:,:,i,j),'salt & pepper', .1);
    end
end