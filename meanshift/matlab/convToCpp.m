function [ code ] = convToCpp( M,type,varname)
%CONVTOCPP Summary of this function goes here
%   Detailed explanation goes here
    
   % code = sprintf('vector< vector<double> > M;\n');
   
   if ( nargin < 2 )
       type = 'double';
   end
   if ( nargin < 3 )
       varname = 'm';
   end
   
   
    dim = size(M,2);
    n = size(M,1);
    
    if (dim == 1)
        code = sprintf('%s %s[%d] = ',type,varname,n); 
    else
        code = sprintf('%s %s[%d][%d] = ',type,varname,n,dim);
    end
    code = [code  sprintf('{ ')];
    for i = 1:n
        if (dim > 1 )
            code = [code  sprintf('{')];
        end
        for d = 1:dim
            code = [code sprintf('%e',M(i,d)) ];
            if  (d < dim)
                    code = [code ', ' ];
            end
        end
        if (dim > 1 )
            code = [code sprintf('}')];
        end
        if ( i < n )
            code = [code sprintf(',\n')];
        end
    end
    code = [code  sprintf('};\n')];
    
    %code = [code sprintf('vector< vector<%1$s> > %2$s (%2$sC,%2$sC + sizeof(%2$sC) / sizeof(%1$s));',type,varname)];
    %code = [code sprintf('vector<%1$s> %2$s (%2$sC,%2$sC + sizeof(%2$sC) / sizeof(%1$s));',type,varname)];
    % myints + sizeof(myints) / sizeof(int) );

end

