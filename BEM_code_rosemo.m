
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                        %%%
%%%                         ROSEMO ENGINEERING                             %%%                                    
%%%       Blade Element Momentum Theory Thrust and Torque Calculator       %%%
%%%                         For Wind Turbine                               %%%
%%%                                                                        %%%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc


Nelements       = input ('Number of Elements: ');
bladeNb         = input ('Number of Blades: ');
bladeradius     = input ('Blade Radius (meters): ');
bladehub        = input ('Hub radius (meters): ');
lamda           = input('Lamda : ');
V0              = input('Free stream wind velocity : ');
omega           = lamda*V0/bladeradius;

r(1) = bladehub;

for i = 2:1:Nelements+1
    
    r(i) =  (bladeradius/Nelements)*(i-1);
    
end


for i = 1:1:Nelements+1
    
    chord(i)          = input (' Chord (meters): ');
    local_lamda(i)    = ((r(i)*omega)/V0);
    local_solidity(i) = ((bladeNb*chord(i))/(2*pi*r(i)));
    beta_deg(i)       = 90 - (((2/3)*atan(1/local_lamda(i)))*180/pi);
    beta_rad(i)       = beta_deg(i)*pi/180;
    twist_deg(i)      = input ('Twist angle (degrees): ');
    twist_rad(i)      = twist_deg(i)*pi/180;
    AoA(i)            = twist_deg(i)-beta_deg(i);
end

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%     INITIAL GUESS    VALUES    %%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
  for i = 1:1:Nelements+1;
     
     cl(i) = input ('Lift coefficient ');
         
  end
  

 
 for i = 1:1:Nelements+1;
     
     a(i)    = ((1 + (4*(cos(beta_rad(i)))^2)/(local_solidity(i)*cl(i)*sin(beta_rad(i))))^-1);
     a_pr(i) = (1-(3*a(i)))/((4*a(i))-1);
     
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %        REPEATED ITERATION          %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 sayac=0;
 
 while sayac<Nelements+1
     
     delta =0.1;
     i=sayac+1;
 
     while delta >0.00001 || a_pr(i) < 0
     
     new_beta_rad(i) = atan(( local_lamda(i)*(1+a_pr(i)))/(1-a(i)));
     new_beta_deg(i) = new_beta_rad(i)*180/pi
         
     AoA(i)          = twist_deg(i)-new_beta_deg(i);
     k               = AoA(i);

     fprintf( 'Angle of attack is %0.4f for the %3d .  blade element\n  ',k,i);
     new_cl(i)       = input ('Input New Lift Coefficient  ');
        
     new_a(i)        = ((1 + (4*(cos(new_beta_rad(i)))^2)/(local_solidity(i)*new_cl(i)*sin(new_beta_rad(i))))^-1);
     new_a_pr(i)     = ((local_solidity(i)*new_cl(i))/(4*local_lamda(i)*cos( new_beta_rad(i))))*(1-new_a(i));
     
     delta_a(i)      = new_a(i)-a(i);
     delta_a_pr(i)   = new_a_pr(i)-a_pr(i);
     
     delta           = abs(delta_a_pr(i))
     
     a(i)            = new_a(i);
     a_pr(i)         = new_a_pr(i);
     
     end
    
     sayac=sayac+1;
 end
   
 for i = 1:1:Nelements+1
     
     K(i)=(local_lamda(i)^3)*a_pr(i)*(1-a(i));
     
 end
 
  for i = 1:1:Nelements
      
     H(i)=((local_lamda(i+1)-local_lamda(i)))/2;
     
  end
  
  curve_area=0;
  
  for i = 1:1:Nelements
      
      integral(i) = (K(i)+K(i+1))*H(i);
      curve_area  = curve_area + integral(i);
      
  end
  
      Cp          = curve_area*(1/8);
  
  