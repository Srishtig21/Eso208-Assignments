function [y] = polynomial (x,B ,n,m)
  P = ones(n,1);
  for i=1:n
    for j=1:n
      if(i~=j)
        P(i,1)=P(i,1).*(x-B(j,1))/(B(i,1)-B(j,1));
      end;
    end;
  end;
  y=0;
  for i=1:n
    y=y+B(i,2)*P(i,1);
  end;
end