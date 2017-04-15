deff('u=f1(y)','u=-sin(y)*sin((%pi*y)/2)');
deff('u=f2(y)','u=sin(y)*sin((3*%pi*y)/2)');
deff('u=f3(x)','u=-sin(x)*sin((%pi*x)/2)');
deff('u=f4(x)','u=sin(x)*sin((3*%pi*x)/2)');
deff('u=f(x,y)','u=-cos(x+y)*sin(x*y)*(2+x^2+y^2)-2*sin(x+y)*cos(x*y)*(x+y)');
deff('z=u(x,y)','z=cos(x+y)*sin(x*y)');

function ks2 = RU(ks1,N1,N2,t,a,b,c,d)
    L1 = b-a;
    L2 = d-c;
    h1 = L1/N1;
    h2 = L2/N2;
    
    k1 = t/h1^2;
    k2 = t/h2^2;
    k3 = 1/(1 + k1 + k2);
    k4 = (1 - k1 - k2);
    ks2 = zeros(N1+1,N2+1);
    
    for i=1:N1+1
        g=(i-1)*h1;
        ks2(i,1) = f3(a+g);
    end

    for j=1:N2+1
        g=(j-1)*h2;
        ks2(1,j) = f1(c+g);
    end
    
    for j=2:N2
        for i=2:N1
            g1 = (i-1)*h1;
            g2 = (j-1)*h2;
            ks2(i,j) = k3 * ( t * ((ks1(i+1,j)+ks2(i-1,j))/h1^2 + (ks1(i,j+1)+ks2(i,j-1))/h2^2 - f(a+g1,c+g2)) + ks1(i,j)*k4);   
        end
    end
    
    for i=1:N1+1
        g=(i-1)*h1;
        ks2(i,N2+1) = f4(a+g);
    end
    
    for j=1:N2+1
        g=(j-1)*h2
        ks2(N1+1,j) = f2(c+g);
    end
endfunction

function ks1 = LD(Un,N1,N2,t,a,b,c,d)
    L1=b-a;
    L2=d-c;
    h1 = L1/N1;
    h2 = L2/N2;
    
    k1 = t/h1^2;
    k2 = t/h2^2;
    k3 = 1/(1 + k1 + k2);
    k4 = (1 - k1 - k2);
    ks1 = zeros(N1+1,N2+1);
    
    for i=1:N1+1
        g=(i-1)*h1;
        ks1(i,N2+1) = f4(a+g);
    end
    
    for j=1:N2+1
        g=(j-1)*h2
        ks1(N1+1,j) = f2(c+g);
    end
    
    for j=N2:-1:2
        for i=N1:-1:2
            g1 = (i-1)*h1;
            g2 = (j-1)*h2;
            ks1(i,j) = k3 * ( t * ((ks1(i+1,j)+Un(i-1,j))/h1^2 + (ks1(i,j+1)+Un(i,j-1))/h2^2 - f(a+g1,c+g2)) + Un(i,j)*k4); 
        end
    end
    
    for i=1:N1+1
        g=(i-1)*h1;
        ks1(i,1) = f3(a+g);
    end

    for j=1:N2+1
        g=(j-1)*h2;
        ks1(1,j) = f1(c+g);
    end
endfunction

a = %pi/2;
b = 3*%pi/2;
c = a;
d = b;
L1 = b-a;
L2 = d-c;
//C = 0;

e = 0.01;


step = [0];
pogStep = [0];
pogUn  = [0];
NN = [0];
tt = [0];
CC = -100:1:100;
for ll = 1:length(CC);
    N1 = 5;
    C = CC(ll);
//    t = tv(ll);
    for l=1:1 
        N1 = 2*N1;
        N2 = N1;
        h = (b-a)/N1;
        
        h1 = L1/N1;
        h2 = L2/N2;
        
        Top = h^2/(sin(%pi*h));
        t = Top;
        n = ceil(N1*log(1/e)/(%pi));
//        printf("n=%d, t=%f, Top=%f\n",n,t,Top);
        
        Un = zeros(N1+1,N2+1);
        Uex = zeros(N1+1,N2+1);
        for i=1:N1+1
            for j=1:N2+1
                g1 = (i-1)*h1;
                g2 = (j-1)*h2;
                Un(i,j) = C;
                Uex(i,j) = u(a+g1,c+g2);
            end
        end
        
        for k=0:1000
            ks1 = LD(Un,N1,N2,t,a,b,c,d);
            ks2 = RU(ks1,N1,N2,t,a,b,c,d);
            
            ks = ks2-ks1;
            Un = ks2;
            pog = max(abs(ks./t));
            pog2 = max(abs(Uex-Un));
                
            if pog <= e  then
//                printf("OK k = %d\n",k);
//                printf("pog=%f, pog2=%f\n",pog,pog2);
                pogStep(l) = pog;
                pogUn(ll) = pog2;
                step(ll) = k;
                NN(ll) = N1;
                tt(ll) = t;
                break;
            end
        end
        x = a:h1:b;
        y = c:h2:d;
    //    plot(x,Uex(:,N2/2)','-b');
    //    plot(x,Un(:,N2/2)','*r');
    //    plot3d1(x,y,Un);
    end
end
