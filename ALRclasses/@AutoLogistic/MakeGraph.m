function G = MakeGraph(r,c,k)
%A static method to produce regular graphs with r*c nodes arranged in r
%rows and c columns, with the number of neighbours of each node
%controlled by k.

%NB: this method isn't very efficient, as it loops through each node in
%turn.  It takes about 3 seconds for a 200x200 graph. Maybe there's a
%better way... or use C++?

q = (1:r*c)';
rn = repmat((1:r)',c,1); %row numbers
cn = repelem((1:c)',r);  %col numbers

    function ix = getidx(i,k)
        %This function gets the indices of nodes j that are within manhattan
        %distance k of node i, with j>i.
        rowdist = abs(rn-rn(i));
        coldist = abs(cn-cn(i));
        ix = find( (rowdist+coldist<=k) & q>i);
    end

%Call the function for each node and build up the pairs of adjacent nodes
%in vectors s and t. Preallocate t and s for speed. What size to
%preallocate to? do a first check on an interior node to see.
%There's a theorem that (sum of degrees of all vertices) = 2*(# edges).
%So get the degree of an interior node and use that to bound the #edges.
intnode = ceil(c/2) * r + ceil(r/2);
intdegree = length(find(abs(rn-rn(intnode)) + abs(cn-cn(intnode)) <= k));
t = zeros(ceil(intdegree*r*c/2),1);
s = zeros(ceil(intdegree*r*c/2),1);
nowrow = 1;
for i = 1:r*c
    tt = getidx(i,k);
    len = length(tt);
    t(nowrow:nowrow+len-1) = tt;
    s(nowrow:nowrow+len-1) = i;
    nowrow = nowrow + len;
end
t(nowrow:end) = [];
s(nowrow:end) = [];
G = graph(s,t);

end
