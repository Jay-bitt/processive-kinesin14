function [M1,m1,h1,state,lengthx,x,x1,x2,x3,x4,x5,x21,x22,x31,x32] = testM(M1,m1,h1,state,N,N1,dn0,lengthx,x,x1,x2,x3,x4,x5,x21,x22,x31,x32,Np)
dt = 1e-6;
kcat = 9*dt; kon = 10*dt; kd = 300*dt; kdiff = 0*dt; k0 = 3700*dt; E = 5; K = 0.2; Fx = 0; r0 = 62.2;
K = K/4.1; 
n = length(m1); n2 = round(N1/dn0 - 1);
if N < N1
    dN = (N-N1)/2; dn = 0;
else
    dN = (N-N1)/2; dn = n2;
end
dN = round(dN);
for i = 1:n
        
    if state(i+dn) == 0  %% combine
        if N1-(m1(i)-dN+round(x/8)) > 0 && N1-(m1(i)+round(x/8)-dN) < N1 && rand() < kon && M1(N1-(m1(i)+round(x/8)-dN)) < 10 
            state(i+dn) = 22; h1(i) = N1-(m1(i)+round(x/8)-dN);
            M1(h1(i)) = 22;
            lengthx(i+dn) = 8*(round(x/8) - x/8);
            if 8*(round(x/8) - x/8) == 4
                if rand() < 0.5
                    lengthx(i+dn) = -4;
                end
            end
            ff = find(state~= 0); nf = length(ff);
            ddx = mean(lengthx(ff)) - Fx/(nf*4.1*K); x = x + ddx; x3 = x3 + ddx; lengthx(ff) = lengthx(ff) - ddx;
        else   %% diffuse
            if m1(i)+1 < N && rand() < kdiff% && M(m1(i)+1) == 0 
                m1(i) = m1(i)+1;
            end
            if m1(i)-1 > 0 && rand() < kdiff% && M(m1(i)-1) == 0
                m1(i) = m1(i)-1;
            end
            
        end
        continue
    elseif state(i+dn) == 22
        if rand() < kd
            state(i+dn) = 32;
        elseif m1(i)+1 < N && rand() < kdiff*exp(-0.5*E)% && M(m1(i)+1) == 0 
            state(i+dn) = 21;
            m1(i) = m1(i)+1;
        else
            ff = find(state ~= 0); nf = length(ff);
            length_left = lengthx; length_left(ff) = length_left(ff) - 8; length_left(i+dn) = length_left(i+dn) + 8;
            meanl = (sum(length_left(ff))/nf - Fx/(nf*4.1*K)); position_left = 8 + meanl;
            length_left(ff) = length_left(ff) - meanl; E0 = 0.5*K*(length_left(ff)*length_left(ff)') + 2*E + Fx*position_left/4.1;
            E1 = 0.5*K*(lengthx(ff)*lengthx(ff)');
            if rand() < k0*exp(-0.5*(E0-E1))
                x = x + position_left; x2 = x2 + position_left; x22 = x22 + position_left;
                state(i+dn) = 21;
                lengthx = length_left;
            end
        end
    elseif state(i+dn) == 21
        if rand() < kd
            state(i+dn) = 31;
        elseif m1(i)-1 > 0 && rand() < kdiff*exp(0.5*E)% && M(m1(i)-1) == 0
            state(i+dn) = 22; 
            m1(i) = m1(i)-1;
        else
            ff = find(state ~= 0); nf = length(ff); 
            length_right = lengthx; length_right(ff) = length_right(ff) + 8; length_right(i+dn) = length_right(i+dn) - 8;
            meanr = (sum(length_right(ff))/nf - Fx/(nf*4.1*K)); position_right = -8 + meanr;
            length_right(ff) = length_right(ff) - meanr; E2 =  0.5*K*(length_right(ff)*length_right(ff)') + Fx*position_right/4.1;
            E1 = 0.5*K*(lengthx(ff)*lengthx(ff)');
            if rand() < k0*exp(-0.5*(E2-E1))
                x = x + position_right; x2 = x2 + position_right;  x21 = x21 + position_right;
                state(i+dn) = 22;
                lengthx = length_right;
            end
        end
    elseif state(i+dn) == 32
        if rand() < kcat
            length_right = lengthx; length_left = lengthx;
            ff = find(state~= 0); nf = length(ff);
            length_left(i+dn) = length_left(i+dn) + 8; ddxl = sum(length_left)/nf - Fx/(nf*4.1*K);
            length_left(ff) = length_left(ff) - ddxl;
            El = 0.5*K*sum(length_left(ff).^2-lengthx(ff).^2) + Fx*ddxl/4.1;
            length_right(i+dn) = length_right(i+dn) - 8; ddxr = sum(length_right)/nf - Fx/(nf*4.1*K);
            length_right(ff) = length_right(ff) - ddxr;
            Er = 0.5*K*sum(length_right(ff).^2-lengthx(ff).^2) + Fx*ddxr/4.1;
            r = r0*exp(-0.5*(El-Er)); pf = r/(1+r); rf = rand();
            if rand() < 1/Np || ((h1(i)+1 == N1 || M1(h1(i)+1) > 10 ) && (h1(i)-1 == 0 || M1(h1(i)-1) > 10))
                state(i+dn) = 0; lengthx(i+dn) = 0; M1(h1(i)) = 0; h1(i) = 0;
                ff = find(state~= 0); nf = length(ff);
                if nf == 0
                    continue
                end
                ddx = mean(lengthx(ff)) - Fx/(nf*4.1*K); x = x + ddx; x4 = x4 + ddx; lengthx(ff) = lengthx(ff) - ddx;
%             elseif M1(h1(i)+1) > 10 && M1(h1(i)-1) > 10
            elseif h1(i)-1 > 0 && M1(h1(i)-1) < 10 && (rf < pf || (h1(i)+1 < N1 && M1(h1(i)+1) > 10) )
                state(i+dn) = 22; lengthx = length_left; x5 = x5 + ddxl; x = x + ddxl;
                M1(h1(i)-1) = 22; M1(h1(i)) = 0; h1(i) = h1(i)-1;
            elseif h1(i)+1 < N1 && M1(h1(i)+1) < 10 && (rf > pf || (h1(i)-1 > 0 && M1(h1(i)-1) > 10) )
                state(i+dn) = 22;  lengthx = length_right; x5 = x5 + ddxr; x = x + ddxr;
                M1(h1(i)+1) = 22; M1(h1(i)) = 0; h1(i) = h1(i)+1;
            end
        elseif m1(i)+1 < N && rand() < kdiff*exp(0.5*E)% && M(m1(i)+1) == 0
            state(i+dn) = 31; 
            m1(i) = m1(i)+1;
        else
            ff = find(state ~= 0); nf = length(ff);
            length_left = lengthx; length_left(ff) = length_left(ff) - 8; length_left(i+dn) = length_left(i+dn) + 8;
            meanl = (sum(length_left(ff))/nf - Fx/(nf*4.1*K)); position_left = 8 + meanl;
            length_left(ff) = length_left(ff) - meanl; E0 = 0.5*K*(length_left(ff)*length_left(ff)') + Fx*position_left/4.1;
            E1 = 0.5*K*(lengthx(ff)*lengthx(ff)');
            if rand() < k0*exp(-0.5*(E0-E1))
                x = x + position_left; x2 = x2 + position_left; x32 = x32 + position_left;
                state(i+dn) = 31;
                lengthx = length_left;
            end    
        end
    elseif state(i+dn) == 31
        if rand() < kcat
            length_right = lengthx; length_left = lengthx;
            ff = find(state~= 0); nf = length(ff);
            length_left(i+dn) = length_left(i+dn) + 8; ddxl = sum(length_left)/nf - Fx/(nf*4.1*K);
            length_left(ff) = length_left(ff) - ddxl;
            El = 0.5*K*sum(length_left(ff).^2-lengthx(ff).^2) + Fx*ddxl/4.1;
            length_right(i+dn) = length_right(i+dn) - 8; ddxr = sum(length_right)/nf - Fx/(nf*4.1*K);
            length_right(ff) = length_right(ff) - ddxr;
            Er = 0.5*K*sum(length_right(ff).^2-lengthx(ff).^2) + Fx*ddxr/4.1;
            r = r0*exp(-0.5*(El-Er)); pf = r/(1+r); rf = rand();
            if rand() < 1/Np || ((h1(i)+1 == N1 || M1(h1(i)+1) > 10 ) && (h1(i)-1 == 0 || M1(h1(i)-1) > 10))
                state(i+dn) = 0; lengthx(i+dn) = 0; M1(h1(i)) = 0; h1(i) = 0;
                ff = find(state~= 0); nf = length(ff);
                if nf == 0
                    continue
                end
                ddx = mean(lengthx(ff)) - Fx/(nf*4.1*K); x = x + ddx; x4 = x4 + ddx; lengthx(ff) = lengthx(ff) - ddx;
%             elseif M1(h1(i)+1) > 10 && M1(h1(i)-1) > 10
            elseif h1(i)-1 > 0 && M1(h1(i)-1) < 10 && (rf < pf || (h1(i)+1 < N1 && M1(h1(i)+1) > 10) )
                state(i+dn) = 21; lengthx = length_left; x5 = x5 + ddxl; x = x + ddxl;
                M1(h1(i)-1) = 21; M1(h1(i)) = 0; h1(i) = h1(i)-1;
            elseif h1(i)+1 < N1 && M1(h1(i)+1) < 10 && (rf > pf || (h1(i)-1 > 0 && M1(h1(i)-1) > 10) )
                state(i+dn) = 21;  lengthx = length_right; x5 = x5 + ddxr; x = x + ddxr;
                M1(h1(i)+1) = 21; M1(h1(i)) = 0; h1(i) = h1(i)+1;
            end
        elseif m1(i)-1 > 0 && rand() < kdiff*exp(-0.5*E)% && M(m1(i)-1) == 0
            state(i+dn) = 32; 
            m1(i) = m1(i)-1;
        else
            ff = find(state ~= 0); nf = length(ff);
            length_right = lengthx; length_right(ff) = length_right(ff) + 8; length_right(i+dn) = length_right(i+dn) - 8;
            meanr = (sum(length_right(ff))/nf - Fx/(nf*4.1*K)); position_right = -8 + meanr;
            length_right(ff) = length_right(ff) - meanr; E2 =  0.5*K*(length_right(ff)*length_right(ff)') + 2*E + Fx*position_right/4.1;
            E1 = 0.5*K*(lengthx(ff)*lengthx(ff)');
            if rand() < k0*exp(-0.5*(E2-E1))
                x = x + position_right; x2 = x2 + position_right; x31 = x31 + position_right;
                state(i+dn) = 32;
                lengthx = length_right;
            end
        end
    end
    ff = find(state~= 0); nf = length(ff);
    if nf == 0 || state(i+dn) == 0
        continue
    end
    if rand() < 0.5
        length_left = lengthx;
        length_left(i+dn) = length_left(i+dn) + 8; ddxl = sum(length_left)/nf - Fx/(nf*4.1*K);
        length_left(ff) = length_left(ff) - ddxl;
        El = 0.5*K*sum(length_left(ff).^2-lengthx(ff).^2) + Fx*ddxl/4.1;
        if rand() < kdiff*exp(-0.5*El) && m1(i)-1 > 0 %&& M(m1(i)-1) == 0
            lengthx = length_left; x1 = x1 + ddxl; x = x + ddxl;
            m1(i) = m1(i)-1;
            continue
        end
        length_right = lengthx;
        length_right(i+dn) = length_right(i+dn) - 8; ddxr = sum(length_right)/nf - Fx/(nf*4.1*K);
        length_right(ff) = length_right(ff) - ddxr;
        Er = 0.5*K*sum(length_right(ff).^2-lengthx(ff).^2) + Fx*ddxr/4.1;
        if rand() < kdiff*exp(-0.5*Er) && m1(i)+1 < N %&& M(m1(i)+1) == 0
            lengthx = length_right; x1 = x1 + ddxr; x = x + ddxr;
            m1(i) = m1(i)+1;
            continue
        end
    else
        length_right = lengthx;
        length_right(i+dn) = length_right(i+dn) - 8; ddxr = sum(length_right)/nf - Fx/(nf*4.1*K);
        length_right(ff) = length_right(ff) - ddxr;
        Er = 0.5*K*sum(length_right(ff).^2-lengthx(ff).^2) + Fx*ddxr/4.1;
        if rand() < kdiff*exp(-0.5*Er) && m1(i)+1 < N% && M(m1(i)+1) == 0
            lengthx = length_right; x1 = x1 + ddxr; x = x + ddxr;
            m1(i) = m1(i)+1;
            continue
        end
        length_left = lengthx;
        length_left(i+dn) = length_left(i+dn) + 8; ddxl = sum(length_left)/nf - Fx/(nf*4.1*K);
        length_left(ff) = length_left(ff) - ddxl;
        El = 0.5*K*sum(length_left(ff).^2-lengthx(ff).^2) + Fx*ddxl/4.1;
        
        if rand() < kdiff*exp(-0.5*El) && m1(i)-1 > 0% && M(m1(i)-1) == 0
            lengthx = length_left; x1 = x1 + ddxl; x = x + ddxl;
            m1(i) = m1(i)-1;
            continue
        end
    end
end
