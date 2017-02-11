%% Algorimo Genetico Binario %%%
clc;
clear all;
close all;
nvdec = 8 ;
nvbin = 4 ;
psize = 100;
ngen  = 100;
probmut = 0.00001; % probabilidad de mutacion  
nelit = 1; %round((0.02*psize) + 1) ; %cantidad de individuos elit 2% o el mejor
tol = 1e-5 ; % toleracia
Sol = zeros(psize,1);
Selec = zeros(psize-nelit,1);
      
%% POBLACION INICIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pop = InitPop2(nvbin, 2*nvdec, psize) ;

%% INICIO DEL PROCESO ITERATIVO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on ;
igen = 0 ; Best = 1; cbest =0; LBest=0 ;

% for igen = 1:ngen 
while Best>tol;

  igen = igen + 1;
  for ipop=1:psize ;  Sol(ipop) = Func2(Pop(:,:,ipop), nvdec, nvbin) ;   end

%% Ordena la Población
  [Sol,ix] = sort(Sol,'descend'); %'ascend') ;
  BestInd =  Pop(:,:,ix(psize)) ; 
  Best = min(Sol) ;
%   Pop = Pop(:,:,ix);
  
%% inicio figuras  
Mean = std(Sol) ; %mean(Sol) ; 

two=2  ;  lim = 9/(2^nvbin-1) ;  x = 0.0 ; y = 0.0 ;    
SignX = BinDec2(BestInd(:,1)) ; SignY = BinDec2(BestInd(:,nvdec+1));    
Ent2X = BinDec2(BestInd(:,2)) ; Ent2Y = BinDec2(BestInd(:,nvdec+2));    
Ent1X = BinDec2(BestInd(:,3)) ; Ent1Y = BinDec2(BestInd(:,nvdec+3));
Dec1X = BinDec2(BestInd(:,4)) ; Dec1Y = BinDec2(BestInd(:,nvdec+4));
Dec2X = BinDec2(BestInd(:,5)) ; Dec2Y = BinDec2(BestInd(:,nvdec+5));
Dec3X = BinDec2(BestInd(:,6)) ; Dec3Y = BinDec2(BestInd(:,nvdec+6));
Dec4X = BinDec2(BestInd(:,7)) ; Dec4Y = BinDec2(BestInd(:,nvdec+7));
Dec5X = BinDec2(BestInd(:,8)) ; Dec5Y = BinDec2(BestInd(:,nvdec+8));    
X = round([Ent2X Ent1X Dec1X Dec2X Dec3X Dec4X Dec5X].*lim) ;
Y = round([Ent2Y Ent1Y Dec1Y Dec2Y Dec3Y Dec4Y Dec5Y].*lim) ;    
for inum=1:nvdec-1; x = x + 10^(2-inum)*X(inum) ; y = y + 10^(2-inum)*Y(inum) ; end    
if SignX<(2^nvbin-1)/2  ;  x = (-1)*x ;  end
if SignY<(2^nvbin-1)/2  ;  y = (-1)*y ;  end
% figure(1) ; plot(x,y,'.') ; grid ;
figure(1) ; plot(igen,Best,'.') ; figure(1) ; plot(igen,Mean,'r+') ;
pause(0.01)
fprintf('--- Best solution: %6.8f ,x: %6.8f , y: %6.8f , Mean: %6.8f \n',Best,x,y,Mean); 
%% fin figuras

%     if (Best==LBest||igen==1)
%       cbest = cbest+1 ; 
%     else
%       cbest = 1 ;
%     end
%     LBest = Best ;

%% Scaling process (pag 320 Chandrupatla)
% fh = max(Sol) ; fl = min(Sol) ;
% C = 0.1*fh - 1.1*fl ;
% D = max([1 (fh+C)]) ;
% Sol = (Sol + C)/D  ;

% Establece las probabilidades  
  Sum = sum(Sol) ;  Prob   = Sol/Sum ;  Prob = cumsum(Prob) ;
%   SP = 2 ; % 1 <= SP <= 2
%   Prob = 2 - SP + 2*(SP-1)*((1:1:psize)' - 1)/(psize-1);
%   Prob = cumsum(Prob)./psize ;
  
%% Elige la población ELITE
  Elit = ix(psize-nelit+1:1:psize) ;
  
%% SELECCION (Ruleta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  I=1 ; pos = psize ;
  while I <= psize-nelit
    if Prob(pos) > rand ; Selec(I) = ix(pos) ; I = I+1 ; end
    if pos>1 ; 
        pos = pos-1 ;
    else
        pos = psize ;
    end
  end
  Selec = [Selec ; Elit] ;
  
  Pop = [Pop(:,:,Selec)] ;
  
%% CRUZAMIENTO (Crossover Simple) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ncros = round(rand*(psize-nelit)*0.5) ; % numero de parejas 
  for icros=1:ncros
    ind1 = 1 + round((psize-nelit-1)*rand) ;
    ind2 = ind1 ;
    while ind2 == ind1 ; 
      ind2 = 1 + round((psize-nelit-1)*rand) ;  
    end   
    IndM = Pop(:,:,ind1) ;
    IndF = Pop(:,:,ind2) ;        
    for inum=1:nvdec*2
        Pop(1:2:nvbin, inum, ind1)  = IndF(1:2:nvbin, inum) ;
        Pop(1:2:nvbin, inum, ind2)  = IndM(1:2:nvbin, inum);
    end    
  end
  
%% MUTACION (Step mutation) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for ind=1:psize-nelit 
    if rand <= probmut
%         fprintf('! \n');
      idec = round(rand*15)+1 ;
      ibin = round(rand*3)+1 ;
      if(Pop(ibin,idec,ind)==0) ; 
        Pop(ibin,idec,ind) = 1 ;
      else
        Pop(ibin,idec,ind) = 0 ;
      end
    end 
  end    
  
end  %igen

hold off
