function [H,WR,WI] = dgeev(A)

% Balance the matrix
% ¶Ô¾ØÕó½øÐÐÆ½ºâ´¦Àí

SFMIN1 = 1.002084180004486D-292;
SFMAX1 = 9.979201547673599D+291; 

SMLNUM = 4.008336720017946D-292;
ULP = 	2.220446049250313D-016;

CS = -9.255963134931783D+061;
SN = -9.255963134931783D+061;

EPS = 2.220446049250313D-016;
% DEGBAL
n = size(A,1);
SCALE = ones(n,1);
NOCONV = 0;
while(1)
    for i = 1:n
        C = norm(A(:,i));
        R = norm(A(i,:));
        CA = abs(max(A(:,i)));
        RA = abs(max(A(i,:)));

        if( C==0 || R ==0)
            continue;
        end

        G = R/2;
        F = 1;
        S = C+R;
        while(1)
            if(C >= G)
                break;
            end
            F = F*2;
            C = C*2;
            CA = CA*2;
            R = R /2;
            G = G/2;
            RA =RA/2;
        end
        G = C/2;

        while(1)
            if(G <R)
                break;
            end
            F = F/2;
            C = C/2;
            G = G/2;
            CA = CA/2;
            R = R*2;
            RA =RA*2;
        end
        if (( C+R) >=  S*0.95)
            continue;
        end
        if( F < 1 && SCALE( i ) < 1 ) 
            if( F*SCALE( i ) <= SFMIN1 )
                continue
            end
        end
        if( F > 1 && SCALE( i ) > 1 ) 
            if( SCALE( i ) >= SFMAX1 / F )
                continue;
            end
        end
        G = 1 /F;
        SCALE(i) = SCALE(i) * F; 
        NOCONV = 1;
        A(i,:) = G * A(i,:);
        A(:,i) = F * A(:,i);
    end
    if( ~NOCONV)
        break;
    end
    NOCONV = 0;
end
% disp(A);

% Reduce to upper Hessenberg form
% DGEHD2
TAU = zeros(n,1);
for i = 1:n-1
    
    %DLARFG
    %Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
    [A(min(i+2,n):end,i),A(i+1,i),TAU(i)] = DLARFG(n-i,A(i+1,i),A(min(i+2,n):end,i));
    
    AII = A(i+1,i);
    A(i+1,i) = 1;    
    
    %DLARF
    % n * n-i
    C = A(:,i+1:n);
    V = A(i+1:n,i);
    
    WORK = C*V;    
    C = C -TAU(i)*WORK*V' ;
    
    %back
    A(:,i+1:n) = C;
    
    %DLARF
    % n-i * n-i
    C = A(i+1:n,i+1:n);
    V = A(i+1:n,i);
    
    WORK = C'*V;   
    C = C -TAU(i)*V*WORK' ;
     
    %back
    A(i+1:n,i+1:n) = C;
    
    A(i+1,i) = AII;
end

% disp(A)
% Compute eigenvalues only

%dimesion < 75*75 small matrix
% ==== Small matrix ====
WR = zeros(n,1);
WI = zeros(n,1);

for j = 1:n-3
    A(j+2,j) = 0;
    A(j+3,j) = 0;
end
if( 1<=n-2)
    A(n,n-2) = 0;
end
ITMAX = 30* max(10,n);
V = zeros(3,1);
H  = A;
DAT1 = 0.75;
DAT2 = -0.4375;
i = n;
L =1;
count =1;
for ITS = 0:ITMAX
%    ITS
%    if(ITS==7)
%        disp(H);
%    end 
   if(i<1)
        return;
    end
    k  = i;
    while k >= 1 +1
        if( abs(H( k, k-1 ) ) <= SMLNUM)
            break;
        end
        TST = abs( H( k-1, k-1 ) ) + abs( H( k, k ) );
        if( TST == 0 ) 
            if( k-2 >= 1 )
                TST = TST + abs( H( k-1, k-2 ) );
            end
            if( k+1 <= n )
                TST = TST + abs( H( k+1, k ) );
            end
        end        
        if ( abs( H( k, k-1 ) ) <= ULP*TST )
               AB = max( abs( H( k, k-1 ) ), abs( H( k-1, k ) ) );
               BA = min( abs( H( k, k-1 ) ), abs( H( k-1, k ) ) );
               AA = max( abs( H( k, k ) ),  abs( H( k-1, k-1 )-H( k, k ) ) );
               BB = min( abs( H( k, k ) ), abs( H( k-1, k-1 )-H( k, k ) ) );
               S = AA + AB;
               if ( BA*( AB / S ) <= max( SMLNUM, ULP*( BB*( AA / S ) ) ) )
                   break;
               end
        end    
        k= k-1;
    end%30    
    L = k;
    if( L > 1)
        H( L, L-1 ) = 0;
    end
    
    if(L >= i-1)
        if( L == i)
            WR( i ) = H( i, i );
            WI(i) = 0;
        elseif( L <= i-1 ) 
            [H( i-1, i-1 ), H( i-1, i ), H( i, i-1 ), H( i, i ), WR( i-1 ), WI( i-1 ), WR( i ), WI( i ),CS,SN] = DLANV2( H( i-1, i-1 ), H( i-1, i ), H( i, i-1 ), H( i, i ), WR( i-1 ), WI( i-1 ), WR( i ), WI( i ), CS,SN);
        end        
        i = L -1;
        continue;
    end
    if (1)
       i1 = L;
       i2 = i;
    end 
    
    if( ITS==10)
        S = abs( H( L+1, L ) ) + abs( H( L+2, L+1 ) );
        H11 = DAT1*S + H( L, L );
        H12 = DAT2*S;
        H21 = S;
        H22 = H11;
    elseif(ITS==20)
        S = abs( H( i, i-1 ) ) + abs( H( i-1, i-2 ) );
        H11 = DAT1*S + H( i, i );
        H12 = DAT2*S;
        H21 = S;
        H22 = H11;
    else
        H11 = H( i-1, i-1 );
        H21 = H( i, i-1 );
        H12 = H( i-1, i );
        H22 = H( i, i );
    end
    S = abs( H11 ) + abs( H12 ) + abs( H21 ) + abs( H22 );
    if( S==0 )
        RT1R = 0;
        RT1I = 0;
        RT2R = 0;
        RT2I = 0;
    else
        H11 = H11 / S;
        H21 = H21 / S;
        H12 = H12 / S;
        H22 = H22 / S;
        TR = ( H11+H22 ) / 2.0;
        DET = ( H11-TR )*( H22-TR ) - H12*H21;
        RTDISC = sqrt( abs( DET ) );
    if( DET >= 0)
        RT1R = TR*S;
        RT2R = RT1R;
        RT1I = RTDISC*S;
        RT2I = -RT1I;
    else
        RT1R = TR + RTDISC;
        RT2R = TR - RTDISC;
        if ( abs( RT1R-H22 ) <= abs( RT2R-H22 ) )
            RT1R = RT1R*S;
            RT2R = RT1R;
        else
            RT2R = RT2R*S;
            RT1R = RT2R;
        end
        RT1I = 0;
        RT2I = 0;
    end
    end
    
    for M = i-2:-1:L
        H21S = H( M+1, M );
        S = abs( H( M, M )-RT2R ) + abs( RT2I ) + abs( H21S );
        H21S = H( M+1, M ) / S;
        V( 1 ) = H21S*H( M, M+1 ) + ( H( M, M )-RT1R )* ( ( H( M, M )-RT2R ) / S ) - RT1I*( RT2I / S );
        V( 2 ) = H21S*( H( M, M )+H( M+1, M+1 )-RT1R-RT2R );
        V( 3 ) = H21S*H( M+2, M+1 );
        S = abs( V( 1 ) ) + abs( V( 2 ) ) + abs( V( 3 ) );
        V( 1 ) = V( 1 ) / S;
        V( 2 ) = V( 2 ) / S;
        V( 3 ) = V( 3 ) / S;
        if M ==L
            break;
        end
        if (abs( H( M, M-1 ) )*( abs( V( 2 ) )+abs( V( 3 ) ) ) <= ULP*abs( V( 1 ) )*( abs( H( M-1, M-1 ) )+abs( H( M,M ) )+abs( H( M+1, M+1 ) ) ) )
            break;
        end
    end
    
    for k = M :i-1
        NR = min(3,i-k+1);
        if(k > M)
            for j = 1:NR
                V(j) = H( k+j-1, k-1 );
            end
        end
        %call DLARFG
        [V(2:end),V(1),T1] = DLARFG(NR,V(1),V(2:end));
        if (k > M)
            H( k, k-1 ) = V( 1 );  
            H( k+1, k-1 ) = 0;
            if ( k <i-1 )
                 H( k+2, k-1 ) =0;
            end
        elseif M > L
            H( k, k-1 ) = H( k, k-1 )*( 1-T1 );
        end
        V2 = V(2);
        T2 = T1* V2;
        if NR ==3
            V3 = V(3);
            T3 = T1*V3;
            for j = k:i2
                SUM = H( k, j ) + V2*H( k+1, j ) + V3*H( k+2, j );
                H( k, j ) = H( k, j ) - SUM*T1;
                H( k+1, j ) = H( k+1, j ) - SUM*T2;
                H( k+2, j ) = H( k+2, j ) - SUM*T3;
            end
            for j = i1:min(k+3,i)
                %count = count +1;
                %if( count==49)
                %    count =count;
                %end
                SUM = H( j, k ) + V2*H( j, k+1 ) + V3*H( j, k+2 );
                H( j, k ) = H( j, k ) - SUM*T1;
                H( j, k+1 ) = H( j, k+1 ) - SUM*T2;
                H( j, k+2 ) = H( j, k+2 ) - SUM*T3;
            end    
        elseif  NR == 2
            for j = k:i2
                SUM = H( k, j ) + V2*H( k+1, j );
                H( k, j ) = H( k, j ) - SUM*T1;
                H( k+1, j ) = H( k+1, j ) - SUM*T2;
            end
            for j = i1:i
                SUM = H( j, k ) + V2*H( j, k+1 );
                H( j, k ) = H( j, k ) - SUM*T1;
                H( j, k+1 ) = H( j, k+1 ) - SUM*T2;
            end
        end
    end
end

end


function  [X,ALPHA,TAU] = DLARFG(N,ALPHA,X)
% SAFMIN = 2.004168360008973D-292;
    if N <=1
        TAU = 0;
        return;
    end
    X = X(:);

    XNORM = norm( X(1 :N-1));
    if (XNORM==0) || abs(XNORM) < 1e-16
    %if (XNORM==0)
        TAU = 0;
    else        
        BETA = -sqrt(ALPHA*ALPHA+ XNORM* XNORM) * sign(ALPHA);
        SAFMIN = 2.004168360008973D-292;
        KNT = 0;       
        if abs(BETA) < SAFMIN && BETA~=0
            RSAFMN = 1.0/ SAFMIN;
            while(1)
                KNT = KNT + 1;
                X(1:N-1) = RSAFMN * X(1:N-1);
                BETA = BETA*RSAFMN;
                ALPHA = ALPHA*RSAFMN;
                if abs(BETA) >= SAFMIN
                    break;
                end
            end
            XNORM = norm( X(1 :N-1));
            BETA = -sign(ALPHA)* sqrt(ALPHA* ALPHA+ XNORM*XNORM);
        end
        assert(BETA~=0);
        TAU = ( BETA-ALPHA ) / BETA;

        X(1:N-1) = 1 / ( ALPHA-BETA ) .* X(1:N-1);
        for j = 1:KNT
            BETA =  BETA * SAFMIN;
        end
        ALPHA = BETA;
    end
    
end


function [ A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN] = DLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
EPS = 2.220446049250313D-016;
if C ==0
    CS =1;
    SN =0;
elseif B==0
     CS = 0;
     SN = 1;
     TEMP = D;
     D = A;
     A = TEMP;
     B = -C;
     C = 0;
elseif ((A-D)==0 && sign(B) ~= sign(C))
    CS =1;
    SN =0;
else
    TEMP = A - D;
     P = 0.5*TEMP;
     BCMAX = max( abs( B ), abs( C ) );
     BCMIS = min( abs( B ), abs( C ) )*sign(B )*sign( C );
     SCALE = max( abs( P ), BCMAX );
     Z = ( P / SCALE )*P + ( BCMAX / SCALE )*BCMIS;

     if Z >= 4*EPS
        Z = P + sign(P)* sqrt( SCALE )*sqrt( Z );
        A = D + Z;
        D = D - ( BCMAX / Z )*BCMIS;
        TAU = sqrt( C*C+Z* Z );
        CS = Z / TAU;
        SN = C / TAU;
        B = B - C;
        C = 0;
     else
         SIGMA = B + C;
         TAU = sqrt( SIGMA*SIGMA+TEMP *TEMP);
         CS = sqrt( 0.5*( 1+abs( SIGMA ) / TAU ) );
         SN = -( P / ( TAU*CS ) )*sign(SIGMA);
         
         AA = A*CS + B*SN;
         BB = -A*SN + B*CS;
         CC = C*CS + D*SN;
         DD  = -C*SN + D*CS;

         A = AA*CS + CC*SN;
         B = BB*CS + DD*SN;
         C = -AA*SN + CC*CS;
         D = -BB*SN + DD*CS;

         TEMP = 0.5*( A+D );
         A = TEMP;
         D = TEMP;

         if C ~= 0
             if B~=0
                 if sign(B) == sign(C) 
                     SAB = sqrt( abs( B ) );
                     SAC = sqrt( abs( C ) );
                     P = sign( SAB*SAC, C );
                     TAU = 1 / sqrt( abs( B+C ) );
                     A = TEMP + P;
                     D = TEMP - P;
                     B = B - C;
                     C = 0;
                     CS1 = SAB*TAU;
                     SN1 = SAC*TAU;
                     TEMP = CS*CS1 - SN*SN1;
                     SN = CS*SN1 + SN*CS1;
                     CS = TEMP;
                 end
             else
              B = -C;
              C = 0;
              TEMP = CS;
              CS = -SN;
              SN = TEMP;
             end
         end
     end
end
             
RT1R = A;
RT2R = D;      

if C ==0
    RT1I = 0;
    RT2I = 0;
else
    RT1I = sqrt( abs( B ) )*sqrt( abs( C ) );
    RT2I = -RT1I;
end

return;
         
         
end

