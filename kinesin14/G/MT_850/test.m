clear
N1 = 1500;
dt = 1e-6; T = 1e7;
dn = 1; Np = 5;
N0 = 125:125:1000; nk = length(N0); num = 0;
nni = 2; vx = zeros(nk,nni); v1 = vx; v2 = vx; v3 = vx; v4 = vx; v21 = vx; v22 = vx; v31 = vx; v32 = vx; v5 = vx;

for ni = 1:nni
    num = 0;
    for N = N0
        n1 = round(N/dn - 1); n2 = round(N1/dn - 1);
        num = num + 1;
        m1 = zeros(1,n1); h1 = m1; state = zeros(1,n1+n2); lengthx = zeros(1,n1+n2);
        M4 = zeros(1,N1); m2 = zeros(1,n2); h2 = m2;
        M3 = zeros(1,N);
        for i = 1:n1
            m1(i) = round(dn*i)+round((rand()-0.5)*dn/2);
        end
        for i = 1:n2
            m2(i) = round(dn*i - dn/2)+round((rand()-0.5)*dn/2);
        end
        position = zeros(1,T); x = 0; x1 = 0; x2 = 0; x3 = 0; x4 = 0; x5 = 0; x22 = 0; x21 = 0; x32 = 0; x31 = 0;
        for i = 1:T
            [M4,m1,h1,state,lengthx,x,x1,x2,x3,x4,x5,x21,x22,x31,x32] = testM(M4,m1,h1,state,N,N1,dn,lengthx,x,x1,x2,x3,x4,x5,x21,x22,x31,x32,Np);
            position(i) = x;
        end
        v = x/(dt*T);
        vx(num,ni) = v; v1(num,ni) = x1; v2(num,ni) = x2; v3(num,ni) = x3; v4(num,ni) = x4; v5(num,ni) = x5;
        v21(num,ni) = x21; v31(num,ni) = x31; v22(num,ni) = x22; v32(num,ni) = x32;
    end
end
now = datetime('now');dateFormat = 'yyyy_mm_dd_HH';dateStr = datestr(now, dateFormat);fileName = ['output_' dateStr '.csv'];
Vx = zeros(nk,6*nni); Vx(:,1:nni) = vx; Vx(:,nni+1:2*nni) = v1; Vx(:,2*nni+1:3*nni) = v2; Vx(:,3*nni+1:4*nni) = v3;
Vx(:,4*nni+1:5*nni) = v4;  Vx(:,5*nni+1:6*nni) = v5;
csvwrite(fileName, Vx);