function [x] = jc(A,b)
	%% Jacobi Method
	%% Solution of x in Ax=b using Jacobi Method
	% * _*Initailize 'A' 'b' & intial guess 'x'*_
	%%
	fprintf('Jacobi Method\n');
	%x=rand(size(b,1),1)*10;
        x=zeros(size(b,1),1);
	n=size(x,1);
	normVal=Inf; 
	[A b]
	x
	%% 
	% * _*Tolerence for method*_

	tol=1e-20; itr=0;
	%% Algorithm: Jacobi Method
	%%
	while normVal>tol
		xold=x;
		
		for i=1:n
		    sigma=0;
		    
		    for j=1:n
		        
		        if j~=i
		            sigma=sigma+A(i,j)*x(j);
		        end
		        
		    end
		    
		    x(i)=(1/A(i,i))*(b(i)-sigma);
		end
		
		itr=itr+1;
		normVal=abs(xold-x);
	end
	%%
	%%
	fprintf('Solution of the system is :\n'); 
	for i=1:n
		fprintf('%f\n',x(i));
	end
	fprintf('\nin %d iterations\n',itr);
end
