function [ Sf, curve ] = jGA1(r, Ps, N, max_Iter, CR, MR)

fun = @jFitnessFunction1; 
funop = str2func(r.evaldir); 
dim = size(r.Y,2);
X   = false(N,dim); 
for i = 1:N
  for d = 1:dim
    if rand() > 0.5
      X(i,d) = true;
    end
  end
end
fit  = zeros(1,N); 
fitG = inf;
for i = 1:N
  fit(i) = fun(r , X(i,:), Ps );
  if funop(fit(i), fitG)
    fitG = fit(i); 
    Xgb  = X(i,:);
  end
end

curve = inf;
t = 1; 
%---Generations start------------------------------------------------
while t <= max_Iter
  Ifit = 1 - fit; 
  prob = Ifit / sum(Ifit);
  X1   = zeros(1,dim);
  X2   = zeros(1,dim);
  z    = 1;
  for i = 1:N
    if rand() < CR
      k1  = jRouletteWheelSelection(prob);
      k2  = jRouletteWheelSelection(prob);
      P1  = X(k1,:);
      P2  = X(k2,:);
      ind = randi([1, dim - 1]);
      X1(z,:) = [P1(1 : ind), P2(ind + 1 : dim)]; 
      X2(z,:) = [P2(1 : ind), P1(ind + 1 : dim)];
      z = z + 1;
    end
  end
  Xnew = [X1; X2];
  Nc   = size(Xnew,1); 
  Fnew = zeros(1,Nc);
  for i = 1:Nc
    for d = 1:dim
      if rand() <= MR
        Xnew(i,d) = 1 - Xnew(i,d);
      end
    end
  end 
  Xnew = logical(Xnew);
  for i = 1:Nc
    Fnew(i) = fun( r , Xnew(i,:), Ps );
    if funop(Fnew(i), fitG)
      Xgb  = Xnew(i,:); 
      fitG = Fnew(i);
    end
  end 
  XX  = [X; Xnew]; 
  FF  = [fit, Fnew]; 
  [FF, idx] = sort(FF,r.optfunc);
  X   = XX(idx(1:N),:); 
  fit = FF(1:N);
  
  curve(t) = fitG; 
  fprintf('\nIteration %d Best (GA)= %f',t,curve(t))
  t = t + 1;
end
Pos   = 1:dim; 
Sf    = Pos(Xgb == 1); 
Nf    = length(Sf);  

function Index = jRouletteWheelSelection(prob)
C = cumsum(prob); 
P = rand();
for i = 1:length(C)
	if C(i) > P
    Index = i;
    break;
  end
end

