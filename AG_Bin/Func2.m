function Sol = Func2(Ind, nvdec, nvbin)
    two=2  ;  lim = 9/(2^nvbin-1) ;  x = 0.0 ; y = 0.0 ;
    
    SignX = BinDec2(Ind(:,1)) ; SignY = BinDec2(Ind(:,nvdec+1));
    
    Ent2X = BinDec2(Ind(:,2)) ; Ent2Y = BinDec2(Ind(:,nvdec+2));    
    Ent1X = BinDec2(Ind(:,3)) ; Ent1Y = BinDec2(Ind(:,nvdec+3));
    Dec1X = BinDec2(Ind(:,4)) ; Dec1Y = BinDec2(Ind(:,nvdec+4));
    Dec2X = BinDec2(Ind(:,5)) ; Dec2Y = BinDec2(Ind(:,nvdec+5));
    Dec3X = BinDec2(Ind(:,6)) ; Dec3Y = BinDec2(Ind(:,nvdec+6));
    Dec4X = BinDec2(Ind(:,7)) ; Dec4Y = BinDec2(Ind(:,nvdec+7));
    Dec5X = BinDec2(Ind(:,8)) ; Dec5Y = BinDec2(Ind(:,nvdec+8));
    
    X = round([Ent2X Ent1X Dec1X Dec2X Dec3X Dec4X Dec5X].*lim) ;
    Y = round([Ent2Y Ent1Y Dec1Y Dec2Y Dec3Y Dec4Y Dec5Y].*lim) ;
    
    for inum=1:nvdec-1
        x = x + 10^(2-inum)*X(inum) ;
        y = y + 10^(2-inum)*Y(inum) ;
    end
    
    if SignX<(2^nvbin-1)/2  ;  x = (-1)*x ;  end
    if SignY<(2^nvbin-1)/2  ;  y = (-1)*y ;  end
    
    Sol = x*x + y*y ; % Paraboloide eliptico
%     Sol = 10*two + (x*x-10*cos(two*pi*x)) + (y*y-10*cos(two*pi*y)); %Rastrigin's function
return ; 
% 	//return 0.01*(x*x + y*y) + pow(sin( x*x + y*y), 2) ;  
% /*** Rastrigin's function ***/
% 	//return 10*two + (x*x-10*cos(two*PI*x)) + (y*y-10*cos(two*PI*y)) ;

% 	//return 0.01*(x*x + y*y) + pow(sin( x*x + y*y), 2) ;  
% /*** Rastrigin's function ***/
% 	//return 10*two + (x*x-10*cos(two*PI*x)) + (y*y-10*cos(two*PI*y)) ;
% }