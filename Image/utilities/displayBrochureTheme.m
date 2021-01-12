function [P, M, select] = displayBrochureTheme()
global select;
select = 0;
P = [];
M = [];
[T, Match] = loadThemes();
n = length(T);
sz = 50;
mg = 20;

grid = 3;
m = ceil(n/grid);
H = sz*m + mg*(m+1);
W = grid*(5*sz + mg) + mg;
data.W = W;
data.H = H;
data.sz = sz;
data.mg = mg;
data.m = m;
data.grid = grid;

I = 0.2*ones(H, W, 3);
for i = 1:m
    is = (i-1)*(sz + mg) + mg;
    for j = 1:grid
        id = (i-1)*grid+j;
        if id > n
            break;
        end
        js = (j-1)*(mg + 5*sz) + mg;
        I(is+1:is+sz,js+1:js+5*sz,:) = createImageFromTheme(T{id}, sz);
    end
end
data.h = figure;
set(data.h, 'Position', [720, 278, W, H]);
set(data.h, 'menubar', 'none');
mh = uimenu(data.h,'Label','File'); 
data.haxes = axes(...    % Axes for plotting the selected plot
                 'Parent', data.h, ...
                 'Units', 'normalized', ...
                 'Position',[0.1 0.1 0.80 0.8]);
axes(data.haxes);
hi = imshow(I);
set(hi, 'ButtonDownFcn',{@chooseTheme_Callback, data}); 
% set(h, 'WindowButtonMotionFcn',{@motion_Callback, data}); 
uiwait(data.h);
if select > 0 && select <= n
    P = T{select};
    M = Match{select};
end


function chooseTheme_Callback(hObject, e, data)
global select
p = get(data.haxes, 'CurrentPoint');
stepH= data.H/data.m;
select = ceil(p(1,2)/stepH);
i = 1;
stepW = data.W/data.grid;
while(~((i-1)*stepW <= p(1,1) && p(1,1) <= i*stepW))
    i = i + 1;
end
select = (select-1)*data.grid + i;
close(data.h);



function I = createImageFromTheme(T, sz)
n = size(T, 1);
I = zeros(sz, sz*n, size(T,2));
for i = 1:n
    for b = 1:size(T,2)
        I(:, (i-1)*sz+1:i*sz, b) = T(i,b);
    end
end

function [T, Match] = loadThemes()
T{1} = [235 235 230; 226 244 33; 30 97 219; 215 24 24; 20 22 21]/255;
Match{1} = [1 2 3 4 5];
T{2} = [70 137 102; 255 240 165; 255 176 59; 182 73 38; 142 40 0]/255;
Match{2} = [2 5 1 4 1];
T{3} = [248 251 251; 18 138 156; 158 240 52; 180 32 212; 20 74 5]/255;
Match{3} = [1 2 3 4 5];
T{4} = [115 96 61; 191 138 73; 242 202 128; 94 90 89; 13 13 13]/255;
Match{4} = [3 2 5 1 5];
T{5} = [185 18 27; 76 27 27; 246 228 151; 252 250 225; 189 141 70]/255;
Match{5} = [4 2 1 5 1];
T{6} = [34 83 120; 22 149 163; 172 240 242; 243 255 223; 235 127 0]/255;
Match{6} = [4 2 5 2 5];
T{7} = [44 62 80; 231 76 60; 236 240 241; 52 152 219; 41 128 185]/255;
Match{7} = [3 4 2 5 2];
T{8} = [133 219 24; 205 232 85; 245 246 212; 167 197 32; 73 63 11]/255;
Match{8} = [2 1 5 4 5];
T{9} = [16 91 99; 255 250 213; 255 211 78; 219 158 54; 189 73 50]/255;
Match{9} = [1 2 3 4 5];
T{10} = [0 67 88; 31 138 112; 190 219 57; 255 225 26; 253 116 0]/255;
Match{10} = [1 2 3 4 5];
T{11} = [232 229 149; 208 168 37; 64 98 124; 38 57 61; 255 250 228]/255;
Match{11} = [1 2 3 4 5];
T{12} = [180 175 145; 120 119 70; 64 65 30; 50 51 29; 192 48 0]/255;
Match{12} = [1 2 3 4 5];
T{13} = [255 248 227; 204 204 159; 51 51 45; 159 180 204; 219 65 5]/255;
Match{13} = [1 2 3 4 5];
T{14} = [46 9 36; 217 0 0; 255 45 0; 255 140 0; 4 117 111]/255;
Match{14} = [1 2 3 4 5];
T{15} = [89 82 65; 184 174 156; 255 255 255; 172 207 204; 138 9 23]/255;
Match{15} = [1 2 3 4 5];
T{16} = [220 53 34; 217 203 158; 55 65 64; 42 44 43; 30 30 32]/255;
Match{16} = [1 2 3 4 5];
T{17} = [0 0 0; 51 51 51; 255 53 139; 1 176 240; 174 238 0]/255;
Match{17} = [1 2 3 4 5];
T{18} = [41 41 41; 91 120 118; 143 158 139; 242 230 182; 65 42 34]/255;
Match{18} = [1 2 3 4 5];
T{19} = [150 202 45; 181 230 85; 237 247 242; 75 181 193; 127 198 188]/255;
Match{19} = [1 2 3 4 5];
T{20} = [28 29 33; 49 53 61; 68 88 120; 146 205 207; 238 239 247]/255;
Match{20} = [1 2 3 4 5];
% T{21} = [90 31 0; 209 87 13; 253 231 149; 71 119 37; 169 204 102]/255;
T{21} = [241 218 26; 181 230 85; 237 247 242; 75 181 193; 127 198 188]/255;
Match{21} = [1 2 3 4 5];