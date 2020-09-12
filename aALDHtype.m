function a = aALDHtype(type)

%ALDH2 activity based on Chan et al, 2020

if type == 1
    a = 1;
elseif type == 2
    a = .015;
elseif type == 3
    a = .6;
elseif type == 4
    a = .325;
elseif type == 5
    a = .36;
elseif type == 6
    a = .125;
elseif type == 7
    a = .23;
else
    fprintf('not a possible isoform')
    





end
