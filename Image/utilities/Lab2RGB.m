function Iout = Lab2RGB(Iin, LUT)
N2 = 256*256;
Iin = round(Iin);
id = Iin(:,1)*N2 + (Iin(:,2)+127)*256 + Iin(:,3) + 128; 
Iout = double(LUT(id,:))/255;
