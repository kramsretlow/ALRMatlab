function ra = assortativity( M )
%ASSORTATITIVITY computes the assortativity coefficient (see Kolaczyk14, p. 66) of
%model M
% **this is quick and dirty***

E = M.Graph.Edges.EndNodes;
lo = M.Coding(1);
hi = M.Coding(2);
ne = size(E,1);

f11 = 0;
f22 = 0;
f12 = 0;

for i = 1:ne
    from = M.Y(E(i,1));
    to = M.Y(E(i,2));
    if from==lo && to==lo
        f11 = f11+1;
    elseif from==hi && to==hi
        f22 = f22+1;
    else
        f12 = f12+1;
    end
end

f11 = f11/ne;
f22 = f22/ne;
f12 = f12/ne/2; %NOTE: need to divide by 2... this isn't in the books...
f1dot = f11 + f12;
f2dot = f22 + f12;

ra = (f11+f22 - f1dot^2 - f2dot^2) / (1 - f1dot^2 - f2dot^2);

end

