function [x]=gs(A,b)

%% Gauss Seidel Method
%% Solution of x in Ax=b using Gauss Seidel Method
% * _*Initailize 'A' 'b' & intial guess 'x'*_
%%
	fprintf('Gauss Seidel Method');
	%x=rand(size(b,1),1)*10;
        x=zeros(size(b,1),1);
	n=size(x,1);
	normVal=Inf; 
	x
	%% 
	% * _*Tolerence for method*_

	tol=1e-10; itr=0;
	%% Algorithm: Gauss Seidel Method
	%%
	while normVal>tol
		x_old=x;
		
		for i=1:n
		    
		    sigma=0;
		    
		    for j=1:i-1
		            sigma=sigma+A(i,j)*x(j);
		    end
		    
		    for j=i+1:n
		            sigma=sigma+A(i,j)*x_old(j);
		    end
		    
		    x(i)=(1/A(i,i))*(b(i)-sigma);
		end
		
		itr=itr+1;
		normVal=norm(x_old-x);
		%normVal
	end
	%%
	fprintf('\nSolution of the system is :\n\n'); 
	for i=1:n
		fprintf('%f\n',x(i));
	end
	fprintf('\nin %d iterations\n',itr);
end
