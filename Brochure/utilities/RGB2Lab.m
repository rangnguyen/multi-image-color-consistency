function Iout = RGB2Lab(Iin, LUT)
Iin = round(Iin*255);
N2 = 256*256;
id = Iin(:,1)*N2 + Iin(:,2)*256 + Iin(:,3) + 1; 
Iout = LUT(id,:);
