function test_broadcast
x       = [-3,-2,-1,0,1,2,3];
y       = [-2,-1,0,1,2];
z       = [-1,0,1];

xv      = reshape(x,[numel(x),1]);
yv      = reshape(y,[1,numel(y)]);
zv      = reshape(z,[1,1,numel(z)]);

[X,Y,Z] = ndgrid(x,y,z);

A       = @(x,y,z) x.^2 + y.^2 + z.^2;

A0      = A(X,Y,Z);
A1      = A(xv,yv,zv);
AI      = A0==A1;

fprintf('TEST: ndgrid == broadcast \t: ')
if sum(AI(:))==numel(AI)
    fprintf('TRUE\n')
else
    fprintf('FALSE\n')
end
    


%% Timing exp-like fun
    function A = expFun(x,y,z,num)
        A = exp(-num*(x.^2 + y.^2 + z.^2));
    end

x   = linspace(-3,3,300);
y   = linspace(-2,2,200);
z   = linspace(-1,1,100);

xv      = reshape(x,[numel(x),1]);
yv      = reshape(y,[1,numel(y)]);
zv      = reshape(z,[1,1,numel(z)]);

[X,Y,Z] = ndgrid(x,y,z);

N_runs  = 100;

% ndgrid
A0      = ones(size(X));
tic
rng(0)
for i=1:N_runs
    num = rand;
    A0  = exp(-num*(X.^2 + Y.^2 + Z.^2)).*A0;
end
t0  = toc;

% broadcast
A0      = ones(size(X));
tic
rng(0)
for i=1:N_runs
    num = rand;
    A0  = exp(-num*(xv.^2 + yv.^2 + zv.^2)).*A0;
end
t1  = toc;

fprintf('TEST: ndgrid exp fun time \t: %4.3fs\n',t0)
fprintf('TEST: broadcast exp fun time \t: %4.3fs\n',t1)


end