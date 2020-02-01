tic
w = 1:1e9;
toc
tic
w = w*4;
toc
tic
for pos = 1:length(w)
    w(pos) = w(pos)*4;
end;
toc