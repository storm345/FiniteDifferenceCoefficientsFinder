%Script to calculate the coefficients an arbitrary stencil for an arbitrary
%order derivative

%Define stencil
M = 2; %Number of datapoints to the left of the point to use
N = 2; %Number of datapoints to the right of the point to use
stencilPts = M+N+1;

%Define dx and order of derivative to find an approximation for
dx = 1; %Change in x between each point in the stencil
derivativeOrder = 4; %Order of derivative to find stencil coefficients for

if stencilPts < derivativeOrder+1 %Matrix system will not be solvable
   error('Insufficient number of stencil points given to approximate required order derivative!'); 
end

J = max(stencilPts,derivativeOrder+1); %Number of rows required to compute unknowns
Mat = coeffsMatrix(M,N,J,dx); %Matrix such that M*a = b where a is vector of stencil coeffs and b is vector of derivative coeffs (beta_n in notes)
b = zeros(J,1);
b(derivativeOrder+1) = 1; %Vector b where elements are 0 except for coeff of derivative want to find

a = linsolve(Mat,b); %Vector containing stencil coefficients
disp("Derivate is approximately: ");
for k=-M:N
    fprintf([num2str(a(k+M+1)),' * f(xi']);
    if k<0
        fprintf('-');
    else
        fprintf('+');
    end
    fprintf([num2str(abs(k)*dx),')']);
    if k==N
        fprintf('\n');
    else
        fprintf('  +  ');
    end
end
disp("Truncation err order: "+(J+1));

%Compute the matrix of taylor series coefficients that once multiplied by a
%vector of the stencil coefficients will give the coefficients (beta_n in
%the notes) of the derivatives which when summed will aproximate our
%derivative
function Mat = coeffsMatrix(M,N,numRows,dx)
    %Function handle to get the coefficient of the nth term (starting from n=0 for the
    %first term) to the taylor series of f(xi + kdx)
    taylorCoeff = @(k,n,dx) ((k*dx)^n) / factorial(n); %From trivial writing out of taylor series

    Mat = zeros(numRows,M+N+1);
    %Populate matrix
    for i=0:numRows-1
        for j=-M:N
            Mat(i+1,j+M+1) = taylorCoeff(j,i,dx);
        end
    end
end