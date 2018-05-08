dx=0.6;
A=[	 1 -1  0  0  0;
	-1  2 -1  0  0;
	 0 -1  2 -1  0;
	 0  0 -1  2 -1;
	 0  0  0 -1  1];
A= A/dx;
b=[7*dx+10 7*dx 7*dx 7*dx 7*dx-31 ]';;
display('Matriz [A b]');
[A b]
x=gs(A,b);
x=jc(A,b);
