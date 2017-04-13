deff('u=f1(y)','u=-sin(y)*sin((%pi*y)/2)');
deff('u=f2(y)','u=sin(y)*sin((3*%pi*y)/2)');
deff('u=f3(x)','u=-sin(x)*sin((%pi*x)/2)');
deff('u=f4(x)','u=sin(x)*sin((3*%pi*x)/2)');
deff('u=f(x,y)','u=-cos(x+y)*sin(x*y)*(2+x^2+y^2)-2*sin(x+y)*cos(x*y)*(x+y)');
deff('u=f1yy(y)','u=sin(y)*sin((%pi*y)/2)*(1+(%pi/2)^2)-%pi*cos(y)*cos((%pi*y)/2)');
deff('u=f2yy(y)','u=-(sin(y)*sin((3*%pi*y)/2)*(1+(3*%pi/2)^2)-3*%pi*cos(y)*cos((3*%pi*y)/2))');
deff('u=f3xx(x)','u=sin(x)*sin((%pi*x)/2)*(1+(%pi/2)^2)-%pi*cos(x)*cos((%pi*x)/2)');
deff('u=f4xx(x)','u=-(sin(x)*sin((3*%pi*x)/2)*(1+(3*%pi/2)^2)-3*%pi*cos(x)*cos((3*%pi*x)/2))');

function d = Dij(Un, N1, N2, t,a,b,c,d,k)
    L1 = b-a;
    L2 = d-c;
    h1 = L1/N1;
    h2 = L2/N2;
    d = zeros(N1,N2);
    
    for i=1:N1
        g = (i-1)*h1;
        d(i,1)= t * (f3xx(a+g) + f(a+g,c));
//        u = t*(f4xx(a+g)+f(a+g,d));
//        printf("u=%f\n",u);
//        d(i,N2)=t*(f4xx(a+g)+f(a+g,d));
        d(i,N2) = t * (f4xx(a+g) + f(a+g,d));
    end
    
    for j=1:N2-1
        g = (j-1)*h2;
        d(1,j) = t * (f1yy(c+g) + f(a,c+g));
        d(N1, j) = t * (f2yy(c+g) + f(b,c+g));
    end
    
    for i=2:N1-1
        for j=2:N2-1
            g1=(i-1)*h1;
            g2=(j-1)*h2;
            d(i,j) = t*((Un(i-1,j) -2*Un(i,j) + Un(i+1,j))/h1^2 +  (Un(i,j-1 ) -2*Un(i,j) + Un(i,j+1))/h2^2 + f(a+g1,c+g2));
            
        end
    end
endfunction

a = %pi/2;
b = 3*%pi/2;
c = a;
d = b;
C = 0;
e = 0.01;

N1 = 10;
N2 = N1;
h = (b-a)/N1;
t = 2*h^2/(sin(%pi*h));
n = ceil(N1*log(1/e)/(2*%pi));
printf("n=%d, t=%f\n",n,t);

Un = zeros(N1,N2);

for i=1:N1
    for j=1:N2
        Un(i,j) = C;
    end
end

d = Dij(UN)
