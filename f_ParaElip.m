function Sol=f_ParaElip(Ind)
    nvars = size(Ind(3:end),2);
    n = Ind(2) ; two=2 ;
    
% 	// Paraboloide eliptico
	Max =  6.5535 ;% // valor máximo del intervalo
	Min = -6.5535 ;% // valor minimo del intervalo	

    IndX = Ind(1:1:n) ; IndY = Ind(n+1:1:nvars) ; 
    
    
	x = Min + (Max-Min)/((two^16)-1)*BinDec(IndX) ; 
	y = Min + (Max-Min)/((two^16)-1)*BinDec(IndY) ; 
    
    Sol = x*x + y*y ;
	return ; % // Paraboloide eliptico
% 	//return 0.01*(x*x + y*y) + pow(sin( x*x + y*y), 2) ;  
% /*** Rastrigin's function ***/
% 	//return 10*two + (x*x-10*cos(two*PI*x)) + (y*y-10*cos(two*PI*y)) ;

