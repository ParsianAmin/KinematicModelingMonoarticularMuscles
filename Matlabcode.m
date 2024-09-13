
clear all
clc

% input

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

syms t LN ;
t0 = 0.0001 ;
u0 = t0 : 0.001 : 1 ;
u = u0*2*pi ;
% u= t0 : pi/48 : 2*pi ; % pi/1536 : pi/12 : pi + pi/1536 ; % 1.57 ;

g_val = 9.8 ;
g = g_val ;

% A Large Number
LN = 10.^ 5 ;

%% generalized coordinates .   .   .   .   .   .   .   .   .   .   .   .


% N = 2 ;
% % 
% q_min = pi/5 ;
% 
% q_max = 4*pi/5 ;
% 
% q_amplitude = (q_max - q_min) / 4;
% 
% % Reduce the amplitude to shrink the range
% 
% q_mean = (q_max + q_min) / 2;
% 
% q_frequency = 1 ;
% 
% for i = 1 : N
%     q(i) = ( q_amplitude * sin(2 * pi * q_frequency * t) + q_mean) ;
% end

% q(1)=100*sin(pi*t);
% q(2)=100*cos(pi*t);




%% bone and muscle structure input

% pm  L  Lc M  I Matrix
pm = ones( 4 , 4 ) ;
N = size (pm , 1);

N = 3 ;
M = 6 ;
% Muscle properties matrix
% Mp  Lo Lto Le Lte R0 M
Mp = ones ( N , 6*N )*0.5;

Z= zeros(N);


%% extracting values from input
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% 
% % L_v = pm(:, 1)' ;
% L_v = [ 0.1 , 0.5 , 0.5 ,0.8 ] ;
% 
% Lc_v = 0.5 *L_v;
% 
% m_v = pm(:, 3)';
% I_v = pm(:, 4)';
% 
% 
% Lo_v = [ 0 , 0 ; 0.075 , 0.05 ; 0.05 , 0.05 ; 0.15 , 0.1 ] ;
% 
% Le_v = [ 0 , 0 ; 0.4 , 0.4 ; 0.35 , 0.35 ; 0.15 , 0.1 ] ;
% 
% R0_v = [ 0 , 0 ; 0.04 , 0.04 ; 0.04 , 0.04 ; 0.1 , 0.05 ] ;


% 2link test
L_v = [ 0 , 0.5 , 0.5 ] ;

Lc_v = 0.5 *L_v;

m_v = pm(:, 3)' ;
I_v = pm(:, 4)' ;

Lo_v = [ 0,0,0,0,0,0; 0,0,0,0,0,0; 0 , 0.1 , 0.2 , 0.2, 0.3 , 0.5 ] ;

Le_v = [ 0,0,0,0,0,0; 0,0,0,0,0,0; 0 , 0.1 , 0.2 , 0.2, 0.5 , 0.3 ] ;

R0_v = [ 0,0,0,0,0,0; 0,0,0,0,0,0; 0 , 0.05 , 0.15 , 0.15 , 0.25 , 0.25 ] ;
% R0_v = 0.75 * R0_v ;
R0_v = R0_v * 0.8 ;

% cs0_v = [ 0,0,0,0,0,0; 0,0,0,0,0,0; 0 , 0.06 , 0.18 , 0.18 , 0.3  , 0.3 ] ;

Lto_v =[ 0,0,0,0,0,0; 0,0,0,0,0,0; 0 , 0.025 , 0.05 , 0.05 , 0.15 , 0.25 ] ; % R0_v .* 0 ; %  Lo_v/10 ; R0_v .* 0 ; %
Lte_v =[ 0,0,0,0,0,0; 0,0,0,0,0,0; 0 , 0.025 , 0.05 , 0.05 , 0.25 , 0.15 ] ; % R0_v .* 0 ; %  Le_v/10 ; R0_v .* 0 ; % 

% Lt_v = R0_v ;


% for i = 1 : N
%     for j = 1 : M
%        if  i == 1
%          
% % %        Lo_v(i,j) = Mp( j , i ).*0.125 ;
% %        Lto_v(i,j) = Lo_v(i,j).*0.01 ; % Mp( j , i+1 ) ;
% %        
% % %        Le_v(i,j) = Mp( j , i+2 ) ;
% %        Lte_v(i,j) = Le_v(i,j).*0.01 ; % Mp( j , i+3 ).*0.5 ;
% %        
% % %        R0_v(i,j) =Mp( j , i+4 ).*0.0625 ;
% %        
%        M_v(i,j) = Mp( j , i+5 ).*2 ;
%             
%        else
%        
% % %        Lo_v(i,j)  =  Mp ( j , 4*(i-1)+1 ) ;           
% %        Lto_v(i,j) = Lo_v(i,j).*0.01 ; %  Mp ( j , 4*(i-1)+2 ).*0.5 ;
% %        
% % %        Le_v(i,j)  = Mp ( j , 4*(i-1)+3 ) ;
% %        Lte_v(i,j) = Le_v(i,j).*0.01 ; %  Mp ( j , 4*(i-1)+4 ).*0.5 ;
% %        
% % %        R0_v(i,j)  = Mp( j , 4*(i-1)+ 5 ).*0.25 ;
% %             
%        M_v(i,j) = Mp( j , 4*(i-1)+ 6 ).*2 ;
%        
%        
%        end  
%     end
% end
% 


% 
% for i = 1 : N
%     for j = 1 : 2
%         
%         qo_v(i,j)=asin( R0_v(i,j) / Lo_v(i,j) );
%         qe_v(i,j)=asin( R0_v(i,j) / Le_v(i,j) );
%         
%         Lmto_v(i,j)= sqrt( (Lo_v(i,j).^2) - (R0_v(i,j).^2) ) ;
%         Lmte_v(i,j)= sqrt( (Le_v(i,j).^2) - (R0_v(i,j).^2) ) ;
%         
%     end
% end

%% symbolic code for inputs

% Loop to make the variables parametric

for i = 1 : N 

    q_t(i) = sym(sprintf('q_t%d', i), 'real');    

    q_s(i) = sym(sprintf('q_s%d', i), 'real');    
    q_sd(i) = sym(sprintf('q_sd%d', i), 'real');
    q_sdd(i) = sym(sprintf('q_sdd%d', i), 'real');
    
end
% 
% for i = 1 : N
%     
%     % for links
%     L(i) = sym(sprintf('L%d', i), 'positive');
%     Lc(i) = sym(sprintf('Lc%d', i), 'positive');
%     m(i) = sym(sprintf('m%d', i), 'positive');
%     I(i) = sym(sprintf('I%d', i), 'positive');
%  
%     for j = 1 : N
% 
%     % for muscle tendon system ( direct properties)
%     
%     % the amount of end and origin position of muscle tendon system on the link
%     Lo(i,j) = sym(sprintf('Lo%d%d', i, j), 'positive');
%     Le(i,j) = sym(sprintf('Le%d%d', i, j), 'positive');
%     
%     % length of origin and end tendon
%     Lto(i,j) = sym(sprintf('Lto%d%d', i, j), 'positive');
%     Lte(i,j) = sym(sprintf('Lte%d%d', i, j), 'positive');
%     
%     R0(i,j) = sym(sprintf('R0%d%d', i, j), 'positive');
%     Cs0(i,j) = sym(sprintf('Cs0%d%d', i, j), 'positive'); 
%     end
% end


L  = L_v  ;
Lc  = Lc_v  ;
I = I_v ;
m = m_v ;
Lo = Lo_v ;
Le = Le_v ;
% Lto = Lto_v ;
% Lte = Lte_v ;
R0 = R0_v ;

% cs0 = cs0_v ;

% the indirect constant inputs
for i = 1 : N
    for j = 1 : M
        
        % the constant angles of origin and end point of MT system
        qo(i,j)  = asin( R0(i,j) / Lo(i,j) );
        qe(i,j)  = asin( R0(i,j) / Le(i,j) );
        
        % the constant length of origin and end part of MT system       
        Lmto(i,j)= sqrt( (Lo(i,j).^2) - (R0(i,j).^2) );
        Lmte(i,j)= sqrt( (Le(i,j).^2) - (R0(i,j).^2) );
        


%         % sum of origin and end tendons and angles
%         Lt(i,j)  = Lto(i,j) + Lte(i,j);
        
    end
end


% % defining  counter functions
% for i = 1 : N
%     
%             CF_1(i)= floor ( 1/i ) ;
%     
%             CF_2(i)= floor( (i+1)/i )- 2.*floor( 1/i ) ;
%             
%     for j = 1 : 2
%         
%             CF_pi(i,j) =  (floor(1/i).* floor( (2.*j + (-1).^(j))/ (2.*j) ))+ ( floor( (2.*(i)+1) / (2.*(i) - (1/(i))) ) - 3.* floor(1/(i))) ;
%             
%             CF_pi_A(i,j)= floor( (2.*j + (-1).^(j))/ (2.*j) );
%             
%             CF_pi_A0(i,j)= floor( (2.*j + (-1).^(j))/ (2.*j) );
%             
%             CF_i1(j)= floor (1/i);
% 
%             CF_j1(j)= floor (1/j);
%             
%             CF_j2(j)= 1 -  floor( 1/j ) ;
%             
% 
%     end 
% end

%% angles calculatuion

qoe = qo + qe ;
Loe = Lo + Le ;
qc_max = pi + qoe ;
qc_m = pi - qoe ;


for i = 1 : N
    for j = 1 : M


          qc(i,j) = qoe(i,j) + ((-1).^(j+1)).*q_t(i) ;
          
          qs(i,j) = pi  - abs( q_t(i) ) ; %+(( -1 ).^ (j+1)).*( q_t(i) );
          
          Lmtoe(i,j) = Lmto(i,j)+Lmte(i,j) ;
          
          Lmtc_m(i,j) = R0(i,j).*qc_max(i,j) ;
          
          Lmt_m(i,j) = Lmtoe(i,j) + Lmtc_m(i,j) ;
          
          Rc_m(i,j) = ( R0(i,j)).*sin(0.5*qc_max(i,j)) / (0.5*qc_max(i,j)) ;
          
          Rmtoe_max(i,j) =  ((0.5.*Lmto(i,j).*cos(qo(i,j)) - Lo_v(i,j)) + (0.5.*Lmte(i,j).*cos(qe(i,j)) - Le_v(i,j) ))/2 ;
          
          Rmt_max(i,j) = (Rc_m(i,j).*Lmtc_m(i,j)) + (Rmtoe_max(i,j).*Lmtoe(i,j)) ;
          Rmt_max(i,j) = Rmt_max(i,j) / Lmt_m(i,j) ;
          Rmt_max(i,j) = L_v(i) - abs (Rmt_max(i,j)) ;
          
          Rmt_min(i,j) = L_v(i) - (( Lo_v(i,j) +  Le_v(i,j) )/2) ;
          
    end
end

% for i = 1 : N
%     for j = 2 : M
%           rmax(i,j) = 0.5.* R0(i,j).*qc_max(i,j) ; % sqrt ( cs0(i,j) / pi ) ;
%           Lcreal(i,j) =  R0(i,j).*qc_max(i,j) ;
%           rmin(i,j) = 2.* ( ( 2.*rmax(i,j) - Lcreal(i,j) ) / ( 4 + qc_max(i,j) ) ) ; 
%           Range(i,j) = R0(i,j)/10 ; % rmax(i,j) - rmin(i,j) ;
%           
%           R = abs(( qc(i,j) / qc_max(i,j) )) .* Range(i,j) ;        
%         
%     end
% end

% R0
% rmax
% rmin
% Loe
% 
% RR =  R0 + rmax ;
% R0 = R0 - R ;
%% MT length rate for both conditions when MT is curve and when MT is straight

% inner angle of MT system in straight mode

for i = 1 : N
    for j = 1 : M
        
     % for when the arch is zero
     
     Lmt_s(i,j)= sqrt( (Lo(i,j).^2) + (Le(i,j).^2) - 2.*Lo(i,j).*Le(i,j).*cos(qs(i,j)) ) ;
     
     Lmtc(i,j) = qc(i,j).*R0(i,j) ;
     
    end
end

%      Lto_v =  Lmto - 0.00000001 ;
%      Lte_v = Lmte - 0.00000001 ;

     Lmt_a = Lmto + Lmtc + Lmte ;
     Lto_a = Lto_v ; % Lto_v/2 + 12.5 * (Lmt_a/100) ;
     Lte_a = Lte_v ; % Lte_v/2 + 12.5 * (Lmt_a/100) ;
     Lt_a  = Lte_a + Lto_a  ;   
    
     Lto_s = Lto_v ; % Lto_v/2 + 12.5 * (Lmt_s/100) ;
     Lte_s = Lte_v ; % Lte_v/2 + 12.5 * (Lmt_s/100) ;
     Lt_s  = Lte_s + Lto_s  ;
     

     Lt = Lto_v + Lte_v ; 
%      
% Lmt_max = Lmto + Lmte + Lcreal ;
% 
% for i = 1 : N
%     for j = 1 : M
%         
%         LPt_a(i,j) = ( Lmt_a(i,j) - ( Lto_v(i,j) + Lte_v(i,j) ) ) / ( Lmt_max(i,j) - ( Lto_v(i,j) + Lte_v(i,j) ) ) ;
%         LPt_s(i,j) = ( Lmt_s(i,j) - ( Lto_v(i,j) + Lte_v(i,j) ) ) / (Lmt_max(i,j) - ( Lto_v(i,j) + Lte_v(i,j) ) ) ;
%         
%         Lto_s(i,j) = Lto_v(i,j).*( 1 + 0.5* LPt_s(i,j) ) ;
%         Lte_s(i,j) = Lte_v(i,j).*( 1 + 0.5* LPt_s(i,j) ) ;
%         
%         Lto_a(i,j) = Lto_v(i,j).*( 1 + 0.5* LPt_a(i,j) ) ;
%         Lte_a(i,j) = Lte_v(i,j).*( 1 + 0.5* LPt_a(i,j) ) ;
% 
%    
%     end
% end

    
     Lt_s  = Lte_s + Lto_s  ;
     Lt_a  = Lte_a + Lto_a  ;   



% the indirect constant variables
% Ltec = sym('Ltec', [N 2]) ;


        Ltoo = Lmto ;
        Ltee= Lmte ;
      
        Lmc = Lmtc ;
        Ltce = Lmtc ;
        Ltco = Lmtc ;
        
        Ltoc =  Lto_a - Ltoo ;
        Ltec =  Lte_a - Ltee ;
        
        Ltoe = Lto_a - ( Ltco + Ltoo) ;
        Lteo = Lte_a - ( Ltce + Ltee) ;
        
        Lmco = Lmtc -  Ltec ;
        Lmce = Lmtc -  Ltoc ;
           

% the geometrical variables and muscle lengths

     Lm_s =  Lmt_s - Lt_s ;

     Lcmto_s = Lto_s + (Lm_s/2) ;
     Lcmte_s = Lte_s + (Lm_s/2) ;
     
     Lm_a = Lmt_a - Lt_a ;  
     Lmo = Lmto - Lto_a ;
     Lme = Lmte - Lte_a ;
     
     Lcmto = Lteo + (Lm_a/2) ;
     Lcmte = Ltoe + (Lm_a/2) ;
             

for i = 1 :N
    for j = 1 : M
        
        % the angle that tendon gets on the curve
        qtoc(i,j) = (Ltoc(i,j))/R0(i,j) ;
        qtec(i,j) = (Ltec(i,j))/R0(i,j) ;
        
    end
end


% 
% % the indirect constant variables
% Ltec = sym('Ltec', [N 2]) ;
% 
% for i = 1 :N
%     for j = 1 : 2
%         
%         Ltoo(i,j) = Lmto(i,j) ;
%         Ltee(i,j) = Lmte(i,j) ;
%         
%         Lmo(i,j) = Lmto(i,j) -  Lto(i,j);
%         Lme(i,j) = Lmte(i,j) -  Lte(i,j);
%         
%         Lmtc(i,j) = qc(i,j).*R0(i,j) ;
%         Lmc(i,j) = Lmtc(i,j) ;
%         
%         Ltec(i,j) = Lmtc(i,j) ;
%         Ltec(i,j) = Lmtc(i,j) ;
%         
%         Ltoe(i,j)= Lto(i,j) - ( Ltec(i,j) + Ltoo(i,j)) ;
%         Lteo(i,j)= Lte(i,j) - ( Ltec(i,j) + Ltee(i,j)) ;
%         
%         Lmco(i,j) = Lmtc(i,j) -  Ltec(i,j);
%         Lmce(i,j) = Lmtc(i,j) -  Ltec(i,j);
%         
% 
%     end
% end
% 
% % the geometrical variables and muscle lengths
% for i = 1 : N
%     for j = 1 : 2
%         
%         
%      Lm_s(i,j) =  Lmt_s(i,j) - Lt(i,j) ;
% 
%      Lcmto_s(i,j) = Lto(i,j) + (Lm_s(i,j)/2) ;
%      
%      Lcmte_s(i,j) = Lte(i,j) + (Lm_s(i,j)/2) ;
%      
%      
%      
%      Lm_a(i,j) = Lmt_a(i,j) - Lt(i,j) ;
%      
%      Lmo(i,j) = Lmto(i,j) - Lto(i,j) ;
%      
%      Lme(i,j) = Lmte(i,j) - Lte(i,j) ;
%      
%      
%      Lcmto1(i,j) = Lteo(i,j) + Lm_a(i,j)/2 ;
% 
%      Lcmte1(i,j) = Ltoe(i,j) + Lm_a(i,j)/2 ;
%              
%         % the angle that tendon gets on the curve
%         qtoc(i,j) = (Ltec(i,j))/R0(i,j) ;
%         qtec(i,j) = (Ltec(i,j))/R0(i,j) ;
%         
%      
%     end
% end


%% FORWARD kinematics and path .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .

syms q [1 N] ;

% Define the minimum and maximum angles for vertical motion
q_min = - pi/4;
q_max = pi/4;

% Calculate the amplitude and mean for the sine function
q_amplitude = (q_max - q_min) / 2;
q_mean = (q_max + q_min) / 2;

% Define the frequency of the sine wave
q_frequency = 1 / (2 * pi) ;

% Calculate the joint angle q over time
q(1)= pi/2 ;
q(2) = q_amplitude * sin(2 * pi * q_frequency * t) + q_mean ;



% q_min2 = - 11*pi/20 ;
% q_max2 = 11*pi/20 ;
q_min2 = - pi/2 ; %  0.000001 ; % - pi + 0.000001 ;
q_max2 =   pi/2 ; % 2*pi - 0.000001 ;
q_mean2 = (q_max2 + q_min2) / 2 ;
q_amplitude2 = (q_max2 - q_min2) / 2 ;
% Define the frequency of the sine wave
% q_frequency = 1 / (2 * pi) ;


q(3) = - q_amplitude2 * sin(2 * pi * q_frequency * t) + q_mean ; % ( pi - 0.000001 )*cos(t) 









% 
% % end effector position
% 
% XE = 0 ; % - 0.01 + 0.01*abs(sin(t)) ; % - 0.05 ; % 0.005 * abs ( sin(2*t) )  ;
% 
% YE = 1 - 0.75*abs(sin(t)) ; % 0.62 +  0.38 *sin(t) ;
% 
% 
% q(1) = pi/2 ; % pi/6 ; % deg2rad(30) ;
% 
% % inverse kinematics put n -1 for 4 links
% for i = 2 : N - 1
%     if i == 2
%         % for the second angle
%     q(i+1) = - q(1) + pi - acos((L(i).^2 + L(i+1).^2 + - XE.^2 + - YE.^2) / (2 .* L(i).* L(i+1)))   ;
%     
%     else
%         % for the first angle
%     q(i-1) = - q(1) + atan2(YE, XE) - atan2((L(i).*sin(q(i))), (L(i-1) + L(i)*cos(q(i))))  ;
% 
%     end
% end
% 
% q(4) =   - ( q(3) + q(2) ) ; % 0.5*pi * ( sin(t) - 1 ); % - q(2) ; % - 0.5*pi + abs( sin( t )) ;
% % q(3)= pi - acos((L(2).^2 + L(3).^2 + - XE.^2 + - YE.^2) / (2 .* L(2).* L(3)))   ;
% % q(2) = atan2(YE, XE) - atan2((L(3).*sin(q(3))), (L(2) + L(3)*cos(q(3))))  ;
% 
% % q(4) = 2*pi - q(3) - atan ( 0.2 )/2 ; % - (pi/2) .* abs(cos(t)) - atan ( 0.2 )/2 ;


for i = 1 : N
    if i == 1
        XE1 = L(i)*cos( q(i) ) ;
        YE1 = L(i)*sin( q(i) ) ;
        
    elseif i == 2
        XE2 = L(i-1)*cos(q(i-1) ) + L(i)*cos(q(i-1)+ q(i)) ;
        YE2 = L(i-1)*sin(q(i-1) ) + L(i)*sin(q(i-1)+ q(i)) ;
        
    elseif i == 3
        XE3 = L(i-2)*cos(q(i-2) ) + L(i-1)*cos(q(i-2)+ q(i-1)) + L(i)*cos(q(i-2)+ q(i-1)+ q(i)) ;
        YE3 = L(i-2)*sin(q(i-2) ) + L(i-1)*sin(q(i-2)+ q(i-1)) + L(i)*sin(q(i-2)+ q(i-1)+ q(i)) ;
%     else
%         XE4 = L(i-3)*cos(q(i-3) ) + L(i-2)*cos(q(i-3)+ q(i-2) ) + L(i-1)*cos(q(i-3)+ q(i-2)+ q(i-1)) + L(i)*cos(q(i-3) + q(i-2)+ q(i-1)+ q(i)) ;
%         YE4 = L(i-3)*sin(q(i-3) ) + L(i-2)*sin(q(i-3)+ q(i-2) ) + L(i-1)*sin(q(i-3)+ q(i-2)+ q(i-1)) + L(i)*sin(q(i-3) + q(i-2)+ q(i-1)+ q(i)) ;  
%     
    end
    
    end

% XE1 = subs( XE1 , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% YE1 = subs( YE1 , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% 
% XE2 = subs( XE2 , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% YE2 = subs( YE2 , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% 
% XE3 = subs( XE3 , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% YE3 = subs( YE3 , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);




XE11 = eval ( subs (XE1 , [t], [u]) );
YE11 = eval ( subs (YE1 , [t], [u]) );

XE22 = eval ( subs (XE2 , [t], [u]) );
YE22 = eval ( subs (YE2 , [t], [u]) );

XE33 = eval ( subs (XE3 , [t], [u]) );
YE33 = eval ( subs (YE3 , [t], [u]) );   

% XE44 = eval ( subs (XE4 , [t], [u]) );
% YE44 = eval ( subs (YE4 , [t], [u]) );   


dXE1= diff (XE1 , t) ;
ddXE1= diff (dXE1 , t) ;

dYE1= diff (YE1 , t) ;
ddYE1= diff (dYE1 , t) ;


dXE2= diff (XE2 , t) ;
ddXE2= diff (dXE2 , t) ;

dYE2= diff (YE2 , t) ;
ddYE2= diff (dYE2 , t) ;

dXE3= diff (XE3 , t) ;
ddXE3= diff (dXE3 , t) ;

dYE3= diff (YE3 , t) ;
ddYE3= diff (dYE3 , t) ;




dXE11 = eval ( subs (dXE1 , [t], [u]) );
dYE11 = eval ( subs (dYE1 , [t], [u]) );

ddXE11 = eval ( subs (ddXE1 , [t], [u]) );
ddYE11 = eval ( subs (ddYE1 , [t], [u]) );


dXE22 = eval ( subs (dXE2 , [t], [u]) );
dYE22 = eval ( subs (dYE2 , [t], [u]) );

ddXE22 = eval ( subs (ddXE2 , [t], [u]) );
ddYE22 = eval ( subs (ddYE2 , [t], [u]) );

dXE33 = eval ( subs (dXE3 , [t], [u]) );
dYE33 = eval ( subs (dYE3 , [t], [u]) );

ddXE33 = eval ( subs (ddXE3 , [t], [u]) );
ddYE33 = eval ( subs (ddYE3 , [t], [u]) );

% 
% q = subs(q , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);

q_d = diff(q , t) ;
% q_d = subs(q_d, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);

q_dd = diff(q_d , t) ;
% q_dd = subs(q_dd, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);

q_val = eval ( subs (q' , [t] , [u] )) ;
dq_val = eval ( subs (q_d' , [t] , [u] )) ;
ddq_val = eval ( subs (q_dd' , [t] , [u] )) ;


 q_val_1 = eval ( subs(q(1), t, u)) ;
 q_val_2 = eval ( subs(q(2), t, u)) ;
 q_val_3 = eval ( subs(q(3), t, u)) ;

%% FORWARD KINEMATICS of Links (bones)

% for links

 for i = 1 : N 
 
      PA(i) = sum(q_t(1:i));

      LX(i) = L(i) .* cos(PA(i));
      LY(i) = L(i) .* sin(PA(i));

 end
 
for i = 1 : N 
    if i == 1
        
        RX(i) = Lc(i) * cos(PA(i));
        RY(i) = Lc(i) * sin(PA(i));
        
    else
        
        RX(i) = sum(LX(1:i-1)) + Lc(i) * cos(PA(i));
        RY(i) = sum(LY(1:i-1)) + Lc(i) * sin(PA(i));
        
    end
end

for i = 1 : N 

    % for links

    dPA(i)  = diff(PA(i) , t);    % Calculate angular velocity dPA_i
    ddPA(i) = diff(dPA(i), t);  % Calculate angular acceleration ddPA_i

    dRX(i)  = diff(RX(i) , t);    % Calculate linear velocity dPX_i
    dRY(i)  = diff(RY(i) , t);    % Calculate linear velocity dPY_i

    ddRX(i) = diff(dRX(i), t);  % Calculate linear acceleration ddPX_i
    ddRY(i) = diff(dRY(i), t);  % Calculate linear acceleration ddPY_i

end

%  Range0=subs(Range, q_t , q) ;
% figure,
% for i = 1 : N
% for j = 1 : M
%         % angles
%         subplot(1, 1, 1 ) ; % 2,1, 3 - i Adjust subplot position (i-1)*2 + 3 - j
%         plot( u , subs(Range0(i,j), t, u)) ;
%         hold on
%         title(['angles vs time ']); % Display iteration indices
%         xlabel('time');
%         ylabel('angles');
%         
% end
% end
%% angles of conditions


for i = 1 : N
    for j = 1 : M
        
        qe_s(i,j)=asin((Lo(i,j)/Lmt_s(i,j)).*sin(qs(i,j)));  
        
%         qo_s(i,j)=asin((Le(i,j)/Lmt_s(i,j)).*sin(qmt_s(i,j)));
%         
    end
end



for i = 1 : N
    for j = 1 : M
            

%             if i == 1
% 
%             qmto(i,j) = CF_j2(j).*pi + ((-1).^(j)).* ( qo(i,j)) ;
%             qro(i,j)  = qmto(i,j) + ((-1).^(j+1)).* ( pi/2 ) ;
%             qrc(i,j)  = (CF_pi_A0(i,j).*pi + ((-1).^(j+1)).* q_t(i)  +  qo(i,j) - qe(i,j) )/2  ;
%             qre(i,j)  = q_t(i) + ((-1).^(j)).* ((pi/2)- qe(i,j)) ;
%             qmte(i,j) = q_t(i) + ((-1).^(j+1)).*qe(i,j) ;           
%             
% %             qmto(i,j) =          + ((-1).^(j  )).*( qo(i,j)+ pi ) ;
% %             qro(i,j)  =          + ((-1).^(j  )).*( qo(i,j) + (pi/2) ) ;
% %             qrc(i,j)  = qro(i,j) + ((-1).^(j+1)).*( qc(i,j).*0.5 ) ;  %((-1).^(j)).*(  pi - qc(i,j) ).* 0.5 
% %             qre(i,j)  = PA(i)    + ((-1).^(j+1)).*( qe(i,j) - (pi/2) ) ;
% %             qmte(i,j) = PA(i)    + ((-1).^(j+1)).*( qe(i,j) ) ;
% 
%             else
%                 
            qmto(i,j) = sum(q_t(1:i-1))  + ((-1).^(j  )).*( qo(i,j)+ pi ) ;
            qro(i,j)  = qmto(i,j) + (-1).^(j+1).*( (pi/2) );
            qrc(i,j)  = qro(i,j) + ((-1).^(j+1)).*( qc(i,j).*0.5 ) ;
            qre(i,j)  = sum(q_t(1:i))   + ((-1).^(j+1)).*( qe(i,j) - (pi/2) ) ;
            qmte(i,j) = qre(i,j)    + ((-1).^(j+1)).*( (pi/2) ) ;
%             
%             qmto(i,j) = sum(q_t(1:i-1))  + ((-1).^(j  )).*( qo(i,j)+ pi ) ;
%             qro(i,j)  = sum(q_t(1:i-1))  + ((-1).^(j  )).*( qo(i,j) + (pi/2) ) ;
%             qrc(i,j)  = qro(i,j) + ((-1).^(j+1)).*( qc(i,j).*0.5 ) ;  %((-1).^(j)).*(  pi - qc(i,j) ).* 0.5 
%             qre(i,j)  = PA(i)    + ((-1).^(j+1)).*( qe(i,j) - (pi/2) ) ;
%             qmte(i,j) = PA(i)    + ((-1).^(j+1)).*( qe(i,j) ) ;
%             
%             end 
    end
end
%% forward kinematics for muscles, condition 0



% third aproach


for i = 1 : N
    for j = 1 : M
%         if i == 1
%             
%             Rmt0_x(i,j) = Le(i,j).*cos( PA(i) ) + Lcmte_s(i,j).*cos( q_t(i) + ((-1).^(j+1)).*( qe_s(i,j)- pi ) );
%             Rmt0_y(i,j) = Le(i,j).*sin( PA(i) ) + Lcmte_s(i,j).*sin( q_t(i) + ((-1).^(j+1)).*( qe_s(i,j)- pi ) );
%             
% %             Pmt0(i,j)   = q_t(i) + ((-1).^(j+1)).*(  pi + qe_s(i,j)) ;
%             
%         else
            Rmt0_x(i,j) = Le(i,j).*cos( PA(i) ) + Lcmte_s(i,j).*cos(  PA(i) + pi + ((-1).^(j+1)).*qe_s(i,j));
            Rmt0_y(i,j) = Le(i,j).*sin( PA(i) ) + Lcmte_s(i,j).*sin(  PA(i) + pi + ((-1).^(j+1)).*qe_s(i,j));
            
%             Pmt0(i,j)   =  pi + ((-1).^(j+1)).*qe_s(i,j) ;
            
%         end
    end
end


% sum(LX(1:i-1))
% Lmt_s
% Lcmt_s
% qmt_s
% qm_s
% qo_s
% Rm_s
% Rmt0_x
% Rmt0_y
% Pmt0
%% forward kinematics for muscles, condition 1

% the kinematics

for i = 1 : N
    for j = 1 : M
        

                 Rmt1_x(i,j) = R0(i,j).*cos( qro(i,j) ) + ( Lcmto(i,j) ).* cos( qmto(i,j) );
                 Rmt1_y(i,j) = R0(i,j).*sin( qro(i,j) ) + ( Lcmto(i,j) ).* sin( qmto(i,j) );
                 
%                  Pmto1(i,j)   = qmto(i,j) ;
                 
    end
end

% sum( LX(1:i-1) ) + 

% Rmt1_x
% Pmt1

%% forward kinematics for muscles, condition 2

% the kinematics of the origin part

for i = 1 : N
    for j = 1 : M
    

                 Rmto2_x(i,j) = R0(i,j).*cos( qro(i,j) ) + ( Lmo(i,j)/2 ).* cos( qmto(i,j) );
                 Rmto2_y(i,j) = R0(i,j).*sin( qro(i,j) ) + ( Lmo(i,j)/2 ).* sin( qmto(i,j) );
                 
%                  Pmto2(i,j)   = qmto(i,j) ;
    end
end

% the kinematics of the curve part

for i = 1 : N
    for j = 1 : M
%     if i == 1
%             
%            qA_2(i,j)   =  qc(i,j) - qtec(i,j) ;
%            Rc_2(i,j)   = ( R0(i,j) / (0.5).* qA_2(i,j) ).*sin(  (0.5).* qA_2(i,j) ) ;
%               
%       
%            Rmtc2_x(i,j) = Rc_2(i,j).* cos( qrc(i,j) + (((-1).^(j)).* qtec(i,j)) ) ;
%            Rmtc2_y(i,j) = Rc_2(i,j).* sin( qrc(i,j) + (((-1).^(j)).* qtec(i,j)) ) ;
%            
% %            Pmtc2(i,j)   = CF_pi(i,j).*pi + ((-1).^(j+1)).*( (0.5).* qA_2(i,j) + qo(i,j)) ;
%            
%         else
            
           qa_2(i,j) =  qc(i,j) - qtec(i,j) ;
           Rc_2(i,j) = ( R0(i,j) / (0.5).* qa_2(i,j) ).*sin(  (0.5).* qa_2(i,j) ) ;
%            
%            qrco(i,j) = qro(i,j) + ((-1).^(j+1)).*( qa_2(i,j).*0.5 ) ;              
             qrco(i,j) = qrc(i,j) + ((-1).^(j)).*( qtec(i,j).*0.5 ) ;
              
           Rmtc2_x(i,j) = Rc_2(i,j).* cos( qrco(i,j) ) ;
           Rmtc2_y(i,j) = Rc_2(i,j).* sin( qrco(i,j) ) ;
          

           
%             qA_2(i,j) =  ((-1).^(j+1)).* q_t(i) + qoe(i,j) - qtec(i,j) ;
%            Rc_2(i,j) = ( R0(i,j) / (0.5).* qA_2(i,j) ).*sin(  (0.5).* qA_2(i,j) ) ;
%         
%            Rmtc2_x(i,j) = Rc_2(i,j).* cos( qrc(i,j) + ((-1).^(j)).* qtec(i,j) ) ;
%            Rmtc2_y(i,j) = Rc_2(i,j).* sin( qrc(i,j) + ((-1).^(j)).* qtec(i,j) ) ;
%                   
%            
           
           
           
%            
%     end
    end
end

% finding the middle of them

for i = 1 : N
    for j = 1 : M
        
        Rmt2_x(i,j)= ( Rmto2_x(i,j) .* ( Lmo(i,j) ) + Rmtc2_x(i,j) .* Lmco(i,j) )/ Lm_a(i,j) ;
        Rmt2_y(i,j)= ( Rmto2_y(i,j) .* ( Lmo(i,j) ) + Rmtc2_y(i,j) .* Lmco(i,j) )/ Lm_a(i,j) ;
        
%       Pmt2(i,j)  = ( Pmto2(i,j) .* (  Lmo(i,j) ) + Pmtc2(i,j) .* ( Lmtc(i,j)- (Lte(i,j) - Lmte(i,j) )) )/ Lm_a(i,j) ;
        
    end
end
        
%   Rmt2_x= simplify(Rmt2_x)      
%   Pmt2 = simplify(Pmt2)

%% forward kinematics for muscles, condition 3

% for the origin part

for i = 1 : N
    for j = 1 : M
        
        Rmto3_x(i,j) = Rmto2_x(i,j) ;
        Rmto3_y(i,j) = Rmto2_y(i,j) ;
%         Pmto3(i,j)   = Pmto2(i,j)   ;
                
    end
end

% for the curve part

for i = 1 : N
    for j = 1 : M
%         if i == 1
%            
%            qa_3(i,j)   =  (0.5).* qc(i,j) ;
%            Rc_3(i,j)   = ( R0(i,j)).*sin( qa_3(i,j) ) / qa_3(i,j) ;
% %          Pmtc3(i,j)   = CF_pi(i,j).*pi + ((-1).^(j+1)).*( (0.5).* qA_2(i,j) + qo(i,j)) ;
%            
%            Rmtc3_x(i,j) = Rc_3(i,j).* cos( qrc(i,j) ) ;
%            Rmtc3_y(i,j) = Rc_3(i,j).* sin( qrc(i,j) ) ;
% %          Pmtc3(i,j)   = qrc(i,j) ;
%            
%         else
%            
           qa_3(i,j)   =  (0.5).* qc(i,j) ;
           Rc_3(i,j)   = ( R0(i,j)).*sin( qa_3(i,j) ) / qa_3(i,j) ;
           
           Rmtc3_x(i,j) = Rc_3(i,j).* cos( qrc(i,j) ) ;
           Rmtc3_y(i,j) = Rc_3(i,j).* sin( qrc(i,j) ) ;
%          Pmtc3(i,j)   = qrc(i,j) ;
          
%         end
    end
end

% for the end part

for i = 1 : N
    for j = 1 : M
%             if i == 1 
%                  
%                  Rmte3_x(i,j) = R0(i,j).*cos ( qre(i,j) ) + ( Lme(i,j)/2 ).*cos ( qmte(i,j) );
%                  Rmte3_y(i,j) = R0(i,j).*sin ( qre(i,j) ) + ( Lme(i,j)/2 ).*sin ( qmte(i,j) );
% %                Pmte3(i,j)   = qmte(i,j) ;
%                  
%             else
%                  
                 Rmte3_x(i,j) = R0(i,j).*cos( qre(i,j) ) + ( Lme(i,j)/2 ).* cos( qmte(i,j) );
                 Rmte3_y(i,j) = R0(i,j).*sin( qre(i,j) ) + ( Lme(i,j)/2 ).* sin( qmte(i,j) );
%                Pmte3(i,j)   = qmte(i,j) ;
                 
%             end
    end
end


% test
% 
% for i = 1 : N
%     for j = 1 : 2
%         
%         Lmf(i,j)= Lmo(i,j)+ Lmtc(i,j) + Lme(i,j);
%     end
% end

% finding the middle of them

for i = 1 : N
    for j = 1 : M
        
        Rmt3_x(i,j) = ( Rmto3_x(i,j).*( Lmo(i,j) ) + Rmtc3_x(i,j).*( Lmtc(i,j) ) + Rmte3_x(i,j).*( Lme(i,j) ) )/ Lm_a(i,j) ;
        Rmt3_y(i,j) = ( Rmto3_y(i,j).*( Lmo(i,j) ) + Rmtc3_y(i,j).*( Lmtc(i,j) ) + Rmte3_y(i,j).*( Lme(i,j) ) )/ Lm_a(i,j) ;
%       Pmt3(i,j)   = ( Pmto3(i,j).*(  Lmo(i,j)  ) + Pmtc3(i,j).*(  Lmtc(i,j)  ) + Pmte3(i,j).*(  Lme(i,j)  ) )/ Lmf(i,j) ;
        
    end
end

%% forward kinematics for muscles, condition 4

% the kinematics of the curve part

for i = 1 : N
    for j = 1 : M
%          if i == 1
%             
%            qA_4(i,j)   =  qmt_s(i,j) + qoe(i,j) - qtoc(i,j) ;
%            Rc_4(i,j)   = ( R0(i,j) / (0.5).* qA_2(i,j) ).*sin(  (0.5).* qA_2(i,j) ) ;
%               
%       
%            Rmtc4_x(i,j) = Rc_2(i,j).* cos( qrc(i,j) + ((-1).^(j+1)).* qtec(i,j) ) ;
%            Rmtc4_y(i,j) = Rc_2(i,j).* sin( qrc(i,j) + ((-1).^(j+1)).* qtec(i,j) ) ;
%            
% %          Pmtc4(i,j)   = CF_pi(i,j).*pi + ((-1).^(j+1)).*( (0.5).* qA_2(i,j) + qo(i,j)) ;
%       
% 
%         else
%             
           qa_4(i,j) =  qc(i,j) - qtoc(i,j) ;
%          qA_4(i,j) =  ((-1).^(j+1)).* q_t(i) + qoe(i,j) - qtoc(i,j) ;
           Rc_4(i,j) = ( R0(i,j) / (0.5).* qa_4(i,j) ).*sin(  (0.5).* qa_4(i,j) ) ;
           
%          qrce(i,j) = qre(i,j) + ((-1).^(j)).*( qa_4(i,j).*0.5 ) ;
           qrce(i,j) = qrc(i,j) + ((-1).^(j+1)).*( qtoc(i,j).*0.5 ) ;

           Rmtc4_x(i,j) = Rc_4(i,j).* cos( qrce(i,j) ) ;
           Rmtc4_y(i,j) = Rc_4(i,j).* sin( qrce(i,j) ) ;
          
           
%     end
    end
end

% the kinematics of the end part

for i = 1 : N
    for j = 1 : M
                    
                 Rmto4_x(i,j) = R0(i,j).*cos ( qre(i,j) ) + ( Lmo(i,j)/2 ).*cos ( qmte(i,j) );
                 Rmto4_y(i,j) = R0(i,j).*sin ( qre(i,j) ) + ( Lmo(i,j)/2 ).*sin ( qmte(i,j) );
%                  Pmto4(i,j)   = qmto(i,j) ;

    end
end

% finding the middle of them

for i = 1 : N
    for j = 1 : M
        
        Rmt4_x(i,j)= ( Rmto4_x(i,j) .* ( Lme(i,j) ) + Rmtc4_x(i,j) .* Lmce(i,j) )/ Lm_a(i,j) ;
        Rmt4_y(i,j)= ( Rmto4_y(i,j) .* ( Lme(i,j) ) + Rmtc4_y(i,j) .* Lmce(i,j) )/ Lm_a(i,j) ;
        
%         Pmt4(i,j)  = ( Pmto4(i,j) .* (  Lme(i,j) ) + Pmtc4(i,j) .* ( Lmtc(i,j)- (Lto(i,j) - Lmto(i,j) )) )/ Lm_a(i,j) ;
        
    end
end
        
%   Rmt4_x= simplify(Rmt2_x)      
%   Pmt4 = simplify(Pmt2)

%% forward kinematics for muscles, condition 5

% the kinematics of the end part

for i = 1 : N
    for j = 1 : M
    
                 Rmt5_x(i,j) = R0(i,j).*cos(  qre(i,j) ) + ( Lcmte(i,j) ).* cos(  qmte(i,j)  );
                 Rmt5_y(i,j) = R0(i,j).*sin(  qre(i,j) ) + ( Lcmte(i,j) ).* sin(  qmte(i,j)  );
                 
%                  Pmte5(i,j)   = qmte(i,j) ;
                 
    end
end

% Rmt1_x
% Pmt1


%% shift functions


% for i = 1 : N
%     for j = 1 : 2
%         
%         if i == 1 
%             
%         SF_0(i,j) = floor ( (  ( (-1).^(j+1)).*(q_t(i)) + qoe(i,j) - CF_j1(j).*pi  )/ LN ) ;
% 
%         else
%             
%         SF_0(i,j) = floor ( (   ( (-1).^(j+1)).*(q_t(i))+ qoe(i,j)   )/ LN ) ;
% 
%         end
%         
%     end
% end

for i = 1 : N
    for j = 1 : M
        
        SF_0(i,j) = floor ( qc(i,j)/ LN ) ;
        
    end
end

for i = 1 : N
    for j = 1 : M
                
        SF_1(i,j) = floor ( (   Lto_a(i,j) - Lmto(i,j)   )/ LN ).*floor( - Lteo(i,j) / LN );
        
        SF_2(i,j) = floor ( (   Lmte(i,j) - Lte_a(i,j)   )/ LN ).*floor( -Ltec(i,j) / LN );
%       SF_2(i,j) = floor ( (   Lmte(i,j) - Lte_a(i,j)   )/ LN ); % .*floor( Lte_a(i,j) - ( Lmte(i,j) + Lmtc(i,j) ) / LN );
         
        SF_3(i,j) = floor ( (   Lto_a(i,j) - Lmto(i,j)   )/ LN ).*floor ( (   Lte_a(i,j) - Lmte(i,j)   )/ LN );
        
        SF_4(i,j) = floor ( (   Lmto(i,j) - Lto_a(i,j)   )/ LN ).*floor( - Ltoc(i,j) / LN );
%       SF_4(i,j) = floor ( (   Lmto(i,j) - Lto_a(i,j)   )/ LN ).*floor( Lto_a(i,j) - ( Lmto(i,j) + Lmtc(i,j) ) / LN );

        SF_5(i,j) = floor ( (   Lte_a(i,j) - Lmte(i,j)   )/ LN ).*floor( - Ltoe(i,j) / LN );
        
    end
end


% for i = 1 : N
%     for j = 1 : 2
%                 
%         SF_1(i,j) = floor ( (   Lto(i,j) - Lmto(i,j)   )/ LN ).*floor( - Lteo(i,j) / LN );
%         
%         SF_2(i,j) = floor ( (   Lmte(i,j) - Lte(i,j)   )/ LN ).*floor( - Ltec(i,j) / LN );
%         
%         SF_3(i,j) = floor ( (   Lto(i,j) - Lmto(i,j)   )/ LN ).*floor ( (   Lte(i,j) - Lmte(i,j)   )/ LN );
%         
%         SF_4(i,j) = floor ( (   Lmto(i,j) - Lto(i,j)   )/ LN ).*floor( - Ltec(i,j) / LN );
%  
%         SF_5(i,j) = floor ( (   Lte(i,j) - Lmte(i,j)   )/ LN ).*floor( - Ltoe(i,j) / LN );
%         
%     end
% end

for i = 1 : N
    for j = 1 : M
        
        SF_0s(i,j) = abs(SF_0(i,j)) ; % floor ( - qc(i,j)/ LN ) + 1 ; % SF_0(i,j).*2  ;
        SF_0a(i,j) = SF_0(i,j) + 1 ;

        SF_11(i,j) = (SF_0a(i,j)) .* (SF_1(i,j)) ;
        
        SF_12(i,j) = (SF_0a(i,j)) .* (SF_2(i,j)) ;
        
        SF_13(i,j) = (SF_0a(i,j)) .* (SF_3(i,j)) ;
        
        SF_14(i,j) = (SF_0a(i,j)) .* (SF_4(i,j)) ;
 
        SF_15(i,j) = (SF_0a(i,j)) .* (SF_5(i,j)) ;
        
    end
end
% %% substituting

% the parts

% Rmt0_x  = subs(Rmt0_x , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Rmt0_y  = subs(Rmt0_y , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% 
% Rmto1_x = subs(Rmto1_x, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Rmto1_y = subs(Rmto1_y, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% 
% Rmt2_x  = subs(Rmt2_x , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Rmt2_y  = subs(Rmt2_y , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% 
% Rmt3_x  = subs(Rmt3_x , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Rmt3_y  = subs(Rmt3_y , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% 
% Rmt4_x  = subs(Rmt4_x , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Rmt4_y  = subs(Rmt4_y , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% 
% Rmte5_x = subs(Rmte5_x, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Rmte5_y = subs(Rmte5_y, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);



% % functions
% SF_0  = subs(SF_0, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]) ;
% SF_11  = subs(SF_11, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]) ;
% SF_12  = subs(SF_12, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]) ;
% SF_13  = subs(SF_13, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]) ;
% SF_14  = subs(SF_14, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]) ;
% SF_15  = subs(SF_15, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]) ;
% SF_Lt = subs(SF_Lt,[ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);

% q = subs( q , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);

% kinematics

%  LX = subs(LX, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
%  LY = subs(LY, [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);


LX  = subs(LX , [q_t], [q]) ;
LY  = subs(LY  , [q_t], [q]) ;


%% speed and velocity of muscles

% Lm_s_v  = subs(Lm_s , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Lm_a_v  = subs(Lm_a , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);
% Lm_s_v  = Lm_s ;
% Lm_a_v  = Lm_a ;

Lmt_s_v = subs(Lmt_s  , [q_t], [q]) ;
Lmt_a_v = subs(Lmt_a  , [q_t], [q]) ;

Lm_s_v  = subs(Lm_s  , [q_t], [q]) ;
Lm_a_v  = subs(Lm_a  , [q_t], [q]) ;

% SF_0_v = subs(SF_0 , [ Lc; L; m; I; Lo; Le; Lto; Lte; R0], [ Lc_v; L_v; m_v; I_v; Lo_v; Le_v; Lto_v; Lte_v; R0_v]);

% SF_0_v = SF_0 ;

SF_0_v  = subs(SF_0  , [q_t], [q]) ;

SF_0s = subs(SF_0s  , [q_t], [q]) ;
SF_0a = subs(SF_0a  , [q_t], [q]) ;

for i = 1 : N
    for j = 1 : M
        
        dLm_s_v(i,j)= diff (Lm_s_v(i,j) , t ) ;
        dLm_a_v(i,j)= diff (Lm_a_v(i,j) , t ) ;
        

        dLmt_s_v(i,j)= diff (Lmt_s_v(i,j) , t ) ;
        dLmt_a_v(i,j)= diff (Lmt_a_v(i,j) , t ) ;       
        
        
        
        
    end
end


for i = 1 : N
    for j = 1 : M
        
%        Lm_v(i,j) = Lm_s_v(i,j) .*(SF_0s(i,j)) + Lm_a_v(i,j) .*(SF_0a(i,j)) ;
%        dLm_v(i,j) = dLm_s_v(i,j) .*(SF_0s(i,j)) + dLm_a_v(i,j) .*(SF_0a(i,j)) ;
%        
       Lmt_v(i,j) = Lmt_s_v(i,j) .*(SF_0s(i,j)) + Lmt_a_v(i,j) .*(SF_0a(i,j)) ;
       dLmt_v(i,j) = dLmt_s_v(i,j) .*(SF_0s(i,j)) + dLmt_a_v(i,j) .*(SF_0a(i,j)) ;
        
       
%        Lt(i,j) = (Lt_s(i,j)).*(SF_0s(i,j)) + (Lt_a(i,j)).*(SF_0a(i,j)) ;
       
       
    end
end

% Lt = subs ( Lt , [q_t], [q]) ;

Lmt_v = subs ( Lmt_v , [q_t], [q]) ;
Lm_v = Lmt_v - Lt ;
dLmt_v = subs ( dLmt_v , [q_t], [q]) ;


%% summations of consitions


% derivitives
% 
% for i = 1 : N
%     for j = 1 : M
%         
%         dRmt0_x(i,j) = diff ( Rmt0_x(i,j) , t ) ;
%         dRmt0_y(i,j) = diff ( Rmt0_y(i,j) , t ) ;
%         
% %         ddRmt0_x(i,j)= diff ( dRmt0_x(i,j), t ) ;
% %         ddRmt0_y(i,j)= diff ( dRmt0_y(i,j), t ) ;
% 
%         
%         
%         dRmt1_x(i,j) = diff ( Rmt1_x(i,j) , t ) ;
%         dRmt1_y(i,j) = diff ( Rmt1_y(i,j) , t ) ;
%         
% %         ddRmt1_x(i,j)= diff ( dRmt1_x(i,j), t  ) ;
% %         ddRmt1_y(i,j)= diff ( dRmt1_y(i,j), t  ) ;
% %              
% %         
%         
%         dRmt2_x(i,j) = diff ( Rmt2_x(i,j) , t ) ;
%         dRmt2_y(i,j) = diff ( Rmt2_y(i,j) , t ) ;
%         
% %         ddRmt2_x(i,j)= diff ( dRmt2_x(i,j), t  ) ;
% %         ddRmt2_y(i,j)= diff ( dRmt2_y(i,j), t  ) ;
% %               
% %         
%         
%         dRmt3_x(i,j) = diff ( Rmt3_x(i,j) , t ) ;
%         dRmt3_y(i,j) = diff ( Rmt3_y(i,j) , t ) ;
%         
% %         ddRmt3_x(i,j)= diff ( dRmt3_x(i,j), t  ) ;
% %         ddRmt3_y(i,j)= diff ( dRmt3_y(i,j), t  ) ;
% %              
% %         
%         
%         dRmt4_x(i,j) = diff ( Rmt4_x(i,j) , t ) ;
%         dRmt4_y(i,j) = diff ( Rmt4_y(i,j) , t ) ;
%         
% %         ddRmt4_x(i,j)= diff ( dRmt4_x(i,j), t  ) ;
% %         ddRmt4_y(i,j)= diff ( dRmt4_y(i,j), t  ) ;
% %              
% %         
% %         
%         dRmt5_x(i,j) = diff ( Rmt5_x(i,j) , t ) ;
%         dRmt5_y(i,j) = diff ( Rmt5_y(i,j) , t ) ;
%         
% %         ddRmt5_x(i,j)= diff ( dRmt5_x(i,j), t  ) ;
% %         ddRmt5_y(i,j)= diff ( dRmt5_y(i,j), t  ) ;
% 
%         
%     end
% end

% using the shift functions
for i = 1 : N
    for j = 1 : M
        
       Rmt0_x(i,j) = (SF_0s(i,j)) .* Rmt0_x(i,j) ;
       Rmt0_y(i,j) = (SF_0s(i,j)) .* Rmt0_y(i,j) ;
       
       Rmt1_x(i,j) = SF_11(i,j) .* Rmt1_x(i,j) ; 
       Rmt1_y(i,j) = SF_11(i,j) .* Rmt1_y(i,j) ; 

       Rmt2_x(i,j)  = SF_12(i,j) .* Rmt2_x(i,j)  ;
       Rmt2_y(i,j)  = SF_12(i,j) .* Rmt2_y(i,j)  ;

       Rmt3_x(i,j)  = SF_13(i,j) .* Rmt3_x(i,j)  ;
       Rmt3_y(i,j)  = SF_13(i,j) .* Rmt3_y(i,j)  ;
       
       Rmt4_x(i,j)  = SF_14(i,j) .* Rmt4_x(i,j)  ;
       Rmt4_y(i,j)  = SF_14(i,j) .* Rmt4_y(i,j)  ;
    
       Rmt5_x(i,j) = SF_15(i,j) .* Rmt5_x(i,j)  ;
       Rmt5_y(i,j) = SF_15(i,j) .* Rmt5_y(i,j)  ; 
       
    end
end

% using the shift functions for derivitives

% 
% for i = 1 : N
%     for j = 1 : M
%         
%        dRmt0_x(i,j) = (SF_0s(i,j)) .* dRmt0_x(i,j) ;
%        dRmt0_y(i,j) = (SF_0s(i,j)) .* dRmt0_y(i,j) ;
% %        
% %        ddRmt0_x(i,j) = (SF_0s(i,j)) .* ddRmt0_x(i,j) ;
% %        ddRmt0_y(i,j) = (SF_0s(i,j)) .* ddRmt0_y(i,j) ;
% %               
% 
%        dRmt1_x(i,j)  = SF_11(i,j) .* dRmt1_x(i,j)  ;
%        dRmt1_y(i,j)  = SF_11(i,j) .* dRmt1_y(i,j)  ;
% 
% %        ddRmt1_x(i,j)  = SF_11(i,j) .* ddRmt1_x(i,j)  ;
% %        ddRmt1_y(i,j)  = SF_11(i,j) .* ddRmt1_y(i,j)  ;
% %        
% %        
%        dRmt2_x(i,j)  = SF_12(i,j) .* dRmt2_x(i,j)  ;
%        dRmt2_y(i,j)  = SF_12(i,j) .* dRmt2_y(i,j)  ;
% % 
% %        ddRmt2_x(i,j)  = SF_12(i,j) .* ddRmt2_x(i,j)  ;
% %        ddRmt2_y(i,j)  = SF_12(i,j) .* ddRmt2_y(i,j)  ;
% %        
% 
%        dRmt3_x(i,j)  = SF_13(i,j) .* dRmt3_x(i,j)  ;
%        dRmt3_y(i,j)  = SF_13(i,j) .* dRmt3_y(i,j)  ;
% 
% %        ddRmt3_x(i,j)  = SF_13(i,j) .* ddRmt3_x(i,j)  ;
% %        ddRmt3_y(i,j)  = SF_13(i,j) .* ddRmt3_y(i,j)  ;
% %        
% 
%        dRmt4_x(i,j)  = SF_14(i,j) .* dRmt4_x(i,j)  ;
%        dRmt4_y(i,j)  = SF_14(i,j) .* dRmt4_y(i,j)  ;
% 
% %        ddRmt4_x(i,j)  = SF_14(i,j) .* ddRmt4_x(i,j)  ;
% %        ddRmt4_y(i,j)  = SF_14(i,j) .* ddRmt4_y(i,j)  ;
% %        
% 
%        dRmt5_x(i,j)  = SF_15(i,j) .* dRmt5_x(i,j)  ;
%        dRmt5_y(i,j)  = SF_15(i,j) .* dRmt5_y(i,j)  ;
% % 
% %        ddRmt5_x(i,j)  = SF_15(i,j) .* ddRmt5_x(i,j)  ;
% %        ddRmt5_y(i,j)  = SF_15(i,j) .* ddRmt5_y(i,j)  ;
% 
%     end
% end
% 




% summation
for i = 1 : N
     for j = 1 : M
         

%          if i == 1
%    Rmt_x(i,j)= Rmt0_x(i,j)+ Rmto1_x(i,j) + Rmt2_x(i,j) + Rmt3_x(i,j) + Rmt4_x(i,j) + Rmte5_x(i,j)  ;
%    Rmt_y(i,j)= Rmt0_y(i,j)+ Rmto1_y(i,j) + Rmt2_y(i,j) + Rmt3_y(i,j) + Rmt4_y(i,j) + Rmte5_y(i,j)  ;
%    
%    
%          else
             
   Rmt_x(i,j)= sum( LX(1:i-1) )  + Rmt0_x(i,j)+ Rmt1_x(i,j) + Rmt2_x(i,j) + Rmt3_x(i,j) + Rmt4_x(i,j) + Rmt5_x(i,j) ;
   Rmt_y(i,j)= sum( LY(1:i-1) )  + Rmt0_y(i,j)+ Rmt1_y(i,j) + Rmt2_y(i,j) + Rmt3_y(i,j) + Rmt4_y(i,j) + Rmt5_y(i,j) ;           
             
             
%          end
     end
end

% latex( Rmt_x(1,1))

 Rmt_x  = subs(Rmt_x  , [q_t], [q]) ;
 Rmt_y  = subs(Rmt_y  , [q_t], [q]) ;

     
% for derivitives


% for i = 1 : N
%      for j = 1 : M
%          
%    dRmt_x(i,j)= dRmt0_x(i,j)+ dRmt1_x(i,j)+ dRmt2_x(i,j)+ dRmt3_x(i,j)+ dRmt4_x(i,j)+ dRmt5_x(i,j)  ;
%    dRmt_y(i,j)= dRmt0_y(i,j)+ dRmt1_y(i,j)+ dRmt2_y(i,j)+ dRmt3_y(i,j)+ dRmt4_y(i,j)+ dRmt5_y(i,j)  ;
% 
% %    ddRmt_x(i,j)= ddRmt0_x(i,j)+ ddRmt1_x(i,j)+ ddRmt2_x(i,j)+ ddRmt3_x(i,j)+ ddRmt4_x(i,j)+ ddRmt5_x(i,j)  ;
% %    ddRmt_y(i,j)= ddRmt0_y(i,j)+ ddRmt1_y(i,j)+ ddRmt2_y(i,j)+ ddRmt3_y(i,j)+ ddRmt4_y(i,j)+ ddRmt5_y(i,j)  ;
% 
%      end
% end
% 
% dRmt_x  = subs(dRmt_x  , [q_t], [q]) ;
% dRmt_y  = subs(dRmt_y  , [q_t], [q]) ;

% ddRmt_x  = subs(ddRmt_x  , [q_t], [q]) ;
% ddRmt_y  = subs(ddRmt_y  , [q_t], [q]) ;
% 

Rc  = eval(subs(Rc_3 , [q_t], [q])) ;

%% plotting and data extracting




% Define joint names
joint_names = {'q1_rz', 'q2_rz', 'q3_rz' };

% Open file for writing
fileID = fopen('angles.mot', 'w');

% Write header
fprintf(fileID, 'Angles\n');
fprintf(fileID, 'version=1\n');
fprintf(fileID, 'nRows=%d\n', length(u));
fprintf(fileID, 'nColumns=%d\n', N ); % Number of joints + time column
fprintf(fileID, 'inDegrees=yes\n\n');
fprintf(fileID, 'Units are S.I. units (second, meters, Newtons, ...)\n');
fprintf(fileID, 'Angles are in degrees.\n\n');
fprintf(fileID, 'endheader\n');

% Write column labels
fprintf(fileID, 'time');
for joint = 1:N -1
    fprintf(fileID, '\t%s', joint_names{joint});
end
fprintf(fileID, '\n');

% Loop over time instances and extract angles
for i = 1:length(u)
    t_val = u(i);
    
    % Evaluate the joint angles at the current time
    q_values = subs(q, t, t_val); % Substitute t with the current time value in q
    
    % Extract q1, q2, q3, q4 and process them
    q1_val = double(q_values(1));
    q2_val = double(q_values(2));
    q3_val = double(q_values(3));
%     q4_val = double(q_values(4));
    
    % Apply specific adjustment to q1_rz, q2_rz, q3_rz, q4_rz
%     q1_val = - 2* q1_val  ; % Adjust q1_rz
    q1_val = - q2_val ; % - (q1_val/2) + (pi/2)  ;
    q2_val = - q3_val ;
%     q3_val = q4_val ;

    % Convert all to degrees
    q_deg =  rad2deg([q1_val, q2_val , q3_val ]);

    % Write time and joint angles to file
    fprintf(fileID, '%f', t_val);
    for joint = 1:N -1
        fprintf(fileID, '\t%f', q_deg(joint));
    end
    fprintf(fileID, '\n');
end

% Close the file
fclose(fileID);




% Specify the full file path
% filePath = 'D:/documents/learning/mechanical engineering/masters Applied design/thesis/Coding/Forward_kinematics_Couple/11p20_R0_Lt.csv';
% 
% filePath = 'D:/documents/learning/mechanical engineering/masters Applied design/thesis/Coding/Forward_kinematics_Couple/position.csv';


% Read the CSV file into a table
dataTable = readtable(filePath);

% Display the first few rows of the table to verify it was read correctly
% disp(dataTable(1:5, :));

% Extract time column
time = dataTable.(dataTable.Properties.VariableNames{1});
% time = dataTable.time;
% Extract force columns
force1 = dataTable.(dataTable.Properties.VariableNames{2});
force2 = dataTable.(dataTable.Properties.VariableNames{3});
force3 = dataTable.(dataTable.Properties.VariableNames{4});
force4 = dataTable.(dataTable.Properties.VariableNames{5});
force5 = dataTable.(dataTable.Properties.VariableNames{6});

X1 = dataTable.(dataTable.Properties.VariableNames{7});
X2 = dataTable.(dataTable.Properties.VariableNames{8});
Y1 = dataTable.(dataTable.Properties.VariableNames{9});
Y2 = dataTable.(dataTable.Properties.VariableNames{10});

% Lm1= dataTable.(dataTable.Properties.VariableNames{11});
% Lm2= dataTable.(dataTable.Properties.VariableNames{12});
% Lm3= dataTable.(dataTable.Properties.VariableNames{13});
% Lm4= dataTable.(dataTable.Properties.VariableNames{14});
% Lm5= dataTable.(dataTable.Properties.VariableNames{15});

RXm1= dataTable.(dataTable.Properties.VariableNames{11});
RXm2= dataTable.(dataTable.Properties.VariableNames{12});
RXm3= dataTable.(dataTable.Properties.VariableNames{13});
RXm4= dataTable.(dataTable.Properties.VariableNames{14});
RXm5= dataTable.(dataTable.Properties.VariableNames{15});


RYm1= dataTable.(dataTable.Properties.VariableNames{16});
RYm2= dataTable.(dataTable.Properties.VariableNames{17});
RYm3= dataTable.(dataTable.Properties.VariableNames{18});
RYm4= dataTable.(dataTable.Properties.VariableNames{19});
RYm5= dataTable.(dataTable.Properties.VariableNames{20});


% 
% RYm1= dataTable.(dataTable.Properties.VariableNames{21});
% RYm2= dataTable.(dataTable.Properties.VariableNames{22});
% RYm3= dataTable.(dataTable.Properties.VariableNames{23});
% RYm4= dataTable.(dataTable.Properties.VariableNames{24});
% RYm5= dataTable.(dataTable.Properties.VariableNames{25});


% Convert the table to a matrix (optional)
dataMatrix = table2array(dataTable);

% Extract individual columns from the matrix (optional)
time = dataMatrix(:, 1);

force1 = dataMatrix(:, 2);
force2 = dataMatrix(:, 3);
force3 = dataMatrix(:, 4);
force4 = dataMatrix(:, 5);
force5 = dataMatrix(:, 6);

X1 = dataMatrix(:, 7);
X2 = dataMatrix(:, 8);
Y1 = dataMatrix(:, 9);
Y2 = dataMatrix(:, 10);


% Lm1 = dataMatrix(:, 11);
% Lm2 = dataMatrix(:, 12);
% Lm3 = dataMatrix(:, 13);
% Lm4 = dataMatrix(:, 14);
% Lm5 = dataMatrix(:, 15);


RXm1 = dataMatrix(:, 11);
RXm2 = dataMatrix(:, 12);
RXm3 = dataMatrix(:, 13);
RXm4 = dataMatrix(:, 14);
RXm5 = dataMatrix(:, 15);



RYm1 = dataMatrix(:, 16);
RYm2 = dataMatrix(:, 17);
RYm3 = dataMatrix(:, 18);
RYm4 = dataMatrix(:, 19);
RYm5 = dataMatrix(:, 20);


% RYm1 = dataMatrix(:, 21);
% RYm2 = dataMatrix(:, 22);
% RYm3 = dataMatrix(:, 23);
% RYm4 = dataMatrix(:, 24);
% RYm5 = dataMatrix(:, 25);



% % Specify the path to your .mot file
% filePath2 = 'D:/documents/learning/mechanical engineering/masters Applied design/thesis/Coding/Forward_kinematics_Couple/LLL.csv';
% 
% % Read the CSV file into a table
% dataTable2 = readtable(filePath2);
% tendon = table2array(dataTable2);
% 
% time3 = tendon(:,2);
% 
% Lt1 = tendon(:,3) ;
% Lt2 = tendon(:,4) ;
% Lt3 = tendon(:,5) ;
% Lt4 = tendon(:,6) ;
% Lt5 = tendon(:,7) ;




% u0 = u0 * 720 ;
% time2 = time/(2*pi) ;
% time4 = time3/(2*pi) ;

% 
% figure,
% % subplot(1, 1, 1);
% % plot (time2 , force1 , '--' , 'LineWidth' , 3) ;
% % hold on
% % plot (time2 , force2 , '--' , 'LineWidth' , 3) ;
% % hold on
% % plot (time2 , force3 , '--' , 'LineWidth' , 3) ;
% % hold on
% % plot (time2 , force4 , '--' , 'LineWidth' , 3) ;
% % hold on
% % plot (time2 , force5 , '--' , 'LineWidth' , 3) ;
% % hold on
% for i = 3 : N
%     for j = 2 : M
% 
%     plot(u0, subs(Lmt_v(i,j), t, u), ':', 'LineWidth' , 2 );
%     title('contraction length');
%     hold on
%     
%     ylabel('contraction length');
%     xlabel('time');
% 
% 
%     end
% end

% Lmt1 = subs(Lmt_v(3,1), t, u)

% 
% 
% figure,
% subplot(1, 1, 1);
% plot (time2 , Lm1 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , Lm2 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , Lm3 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , Lm4 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , Lm5 , '--' , 'LineWidth' , 3) ;
% hold on
% for i = 3 : N
%     for j = 2 : M
% 
%     plot(u0, subs(Lm_v(i,j), t, u), ':', 'LineWidth' , 2 );
%     title('contraction length');
%     hold on
%     
%     ylabel('contraction length');
%     xlabel('time');
% 
% 
%     end
% end



% 
% figure,
% subplot(1, 2, 1);
% plot (time2 , -X1 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , -X2 , '--' , 'LineWidth' , 3) ;
% hold on
% plot(u0,  XE22 , ':', 'LineWidth' , 2 );
% hold on
% plot(u0,  XE33 , ':', 'LineWidth' , 2 );
% hold on
%     title('x position vs time');
%     ylabel('X position');
%     xlabel('time');
% subplot(1, 2, 2);
% plot (time2 , Y1 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , Y2 , '--' , 'LineWidth' , 3) ;
% hold on
% plot(u0,  YE22 , ':', 'LineWidth' , 2 );
% hold on
% plot(u0,  YE33 , ':', 'LineWidth' , 2 );
% hold on
%     title('Y position vs time');
%     ylabel('Y position');
%     xlabel('time');



% 
% 
% 
% figure,
% 
% subplot(1, 2, 1);
% plot (time2 , -X1 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , -RXm1 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , -RXm2 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , -RXm3 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , -RXm4 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , -RXm5 , '--' , 'LineWidth' , 3) ;
% hold on
% 
% plot(u0,  XE22 , '--', 'LineWidth' , 3 );
% hold on
% for i = 3 % 2 : N
%     for j = 2 : M
% 
%     plot(u0, subs(Rmt_x(i,j), t, u) , ':', 'Linewidth' , 2 );
%     hold on
% 
%     end
% end
%  title('horizontal position');
%  ylabel('horizontal position');
%  xlabel('time');
% 
%  
%  subplot(1, 2, 2);
% plot(u0,  XE22 , ':', 'LineWidth' , 2 );
% hold on
% plot (time2 , Y1 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , RYm1 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , RYm2 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , RYm3 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , RYm4 , '--' , 'LineWidth' , 3) ;
% hold on
% plot (time2 , RYm5 , '--' , 'LineWidth' , 3) ;
% hold on
% 
% plot(u0,  YE22 , '--', 'LineWidth' , 3 );
% hold on
% for i = 3 % : N
%     for j = 2 : M
% 
%     plot(u0, subs(Rmt_y(i,j), t, u) , ':', 'Linewidth' , 2  );
%     hold on
% 
%     end
% end
% title('vertical position');
% ylabel('vertical position');
% xlabel('time');
% 
% 
% 
% figure,
% 
% for i = 3 % : N
%     for j = 3 : 4
% 
%     plot(u0, subs(Rmt_y(i,j), t, u) , ':', 'Linewidth' , 2  );
%     hold on
% 
% 
%     end
% end
%     plot(u0,  YE22 , 'LineWidth' , 3 );
%     hold on
%     plot(u0, subs(Rmt_min(3,3), t, u) ,  'Linewidth' , 3  );
%     hold on
%     
%     plot(u0, subs(Rmt_max(3,3), t, u) ,  'Linewidth' , 3  );
%     hold on
%     
% title('vertical position');
% ylabel('vertical position');
% xlabel('time');
% 
% 
% 
% 
% 
% figure,
% for i = 2 : N
%     for j = 1 : M
% 
% 
%     subplot(2, 2, 3);
%     plot(u, subs(dRmt_x(i,j), t, u) , ':', 'Linewidth' , 2 );
%     title('horizontal velocity');
%     hold on
%     plot( u, dXE33 , 'g' );
% 
%     ylabel('horizontal velocity');
%     xlabel('time');
% 
%   
%     subplot(2, 2, 4);
%     plot(u, subs(dRmt_y(i,j), t, u)  , ':', 'Linewidth' , 2 );
%     title('vertical velocity');
%     hold on
%     plot( u, dYE33 , 'g' );
%     
%     ylabel('vertical velocity');
%     xlabel('time');
% 
% 
%     end
% end

