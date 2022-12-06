%% Algorimo Genetico Binario %%%
clc; clear all; close all;

fprintf('--------------------- AG Bin -----------------------\n') 
fprintf(' Genetic algorithm code for heuristic optimization  \n')
fprintf(' of generic objective function with BINARY variables \n')
fprintf(' Authors: Dr. Ing. J.L. Mroginski \n')
fprintf(' Version: 3.0 (2021)\n')
fprintf('----------------------------------------------------------------------\n') 

%% DATOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
funct ='f_ParaElip' ; % función a optimizar 
ngen = 2 ; % nro. de genes (variables de la función)
ncro = 2^4 ; %nro. de cromosomas (cantidad de binarios para escribir un gen)
nvars = ncro*ngen ; 
psize = 1000; % tamaño de la poblacion 
ngener  = 100; %  % maximo numero de generaciones (iteraciones) 
probmut = 2; %  probabilidad de mutacion maxima en %
% porpas = 70; % 70 Porcentaje de individuos que pasan
% npass = 40; % numero de individuos que pasan
nelit = 2 ; % cantidad de individuos elit (en %) o solo el mejor
tol = 1e-3 ; % tolerancia de convergencia

Sol = zeros(psize,1); 
nelit = ceil(nelit/100*psize);
Selec = zeros(psize-nelit,1);

%% POBLACION INICIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pop = InitPop(nvars, psize) ;

%% INICIO DEL PROCESO ITERATIVO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on ;
igen = 0 ; Best = 1; cbest =0; LBest=0 ;
 
while (Best>tol||igen<ngener)
    igen = igen + 1;
    for ipop=1:psize 
        Sol(ipop) = feval(funct,[ngen,ncro,Pop(ipop,:)]) ;   
    end

%% Ordena la Poblacion
    [Sol,ix] = sort(Sol,'descend'); 
    BestInd =  Pop(ix(psize),:) ; 
  
%% inicio figuras  
    Best = min(Sol) ;
    Mean = mean(Sol) ;
    fprintf('--- Best solution:  %6.8f , Mean:  %6.8f \n',Best,Mean);  
    figure(1) ; plot(igen,Mean+std(Sol),'bo') ;
    figure(1) ; plot(igen,Mean-std(Sol),'bo') ;
    figure(1) ; plot(igen,Mean,'r+') ;
% %   figure(1) ; plot(igen,abs(Mean-Best),'go') ;
% 	Max =  6.5535         ;	Min = -6.5535 ;
%   IndX = BestInd(1:1:nvars/2) ; IndY = BestInd(nvars/2+1:1:nvars) ;
% 	x = Min + (Max-Min)/((2^16)-1)*BinDec(IndX) ; 
% 	y = Min + (Max-Min)/((2^16)-1)*BinDec(IndY) ;  
%   figure(1) ; plot(x,y,'.') ; grid ;  
%   pause(0.01)
%% fin figuras

%% Establece las probabilidades  
    Sum    = sum(Sol) ;   Prob   = Sol/Sum ;  Prob   = cumsum(Prob) ;

%% Elige la poblacion ELITE
    Elit = Pop(ix(psize-nelit+1:1:psize),:) ;
  
%% SELECCION (Ruleta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pos = psize;
    I=1;
    while I <= psize-nelit
        if Prob(pos) > rand
            Selec(I) = ix(pos) ;
            I = I+1 ;
        end
        if pos>1
            pos = pos-1 ;
        else
            pos = psize ;
        end
    end

    Pop = [Pop(Selec,:) ; Elit] ;
  
%% CRUZAMIENTO  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ncros = round(rand*(psize-nelit)*0.5) ; % numero de parejas 
    for icros=1:ncros
        ind1 = 1 + round((psize-nelit-1)*rand) ;
        ind2= ind1 ;
        while ind2 == ind1 
            ind2 = 1 + round((psize-nelit-1)*rand) ;  
        end   
        IndM = Pop(ind1,:) ;
        IndF = Pop(ind2,:) ; 
        
        if mod(igen,2)==0 % generacion par (mascara de filtrado)
            mask = round(rand(1,nvars)) ; 
            Pop(ind1,find(mask==1)) = IndF(find(mask==1));
            Pop(ind2,find(mask==0)) = IndM(find(mask==0)); 
        else % generacion impar (Crossover Simple)
            Pop(ind1,1:2:nvars)  = IndF(1:2:nvars);
            Pop(ind2,1:2:nvars)  = IndM(1:2:nvars);
        end
    end
  
%%% MUTACION (Step mutation) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ind=1:psize-nelit 
        if rand <= probmut/100 
            ivar = ceil(rand*nvars) ;
            if(Pop(ind,ivar)==0) 
                Pop(ind,ivar) = 1 ;
            else
                Pop(ind,ivar) = 0 ;
            end
        end 
    end  
  
end  %igen 

hold off

fprintf('-------------------------------- THE END -----------------------------\n') 



% %%% Selecciona los individuos ELIT
%   MaxSol = max(Sol) ;
%   Sol    = Sol-MaxSol-epsilon  ;
%   [Sol,ix] = sort(Sol,'descend') ; %'ascend'
%   Elit = [];
%   ieli = PopSize;
%   for inel = 1:Nelit
%     Elit = [Elit ; NewPob(ix(ieli),:)] ;
%     ieli = ieli-1;
%   end
% 
% %%% Establece las probabilidades  
%   Sum    = sum(Sol); %- MinSol*PopSize
%   Prob   = Sol/Sum ;
%   Prob   = cumsum(Prob) ;
% 
% %%% SELECCION (Ruleta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   tiros = sort(rand(1, PopSize)) ;
%   Selec=[];
%   pos = 1;
%   I=1;
%   while I <= PopSize-Nelit
%     if tiros(I) < Prob(pos) 
%       Selec = [Selec; ix(pos)] ;
%       I=I+1 ;
%       pos = 1;
%     else
%       pos = pos+1;
%     end
%   end
% 
% 
% %%% CRUZAMIENTO (Crossover Simple) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %NewPob = InitPob(Selec,:) ;
%   NewPob = NewPob(Selec,:) ;
%   ncross = round(rand*(PopSize-Nelit)*0.5) ; % numero de parejas 
%   for icros=1:ncross
%     ind(1) = 1 + round((PopSize-Nelit-1)*rand) ;
%     ind(2)= ind(1) ;
%     while ind(2) == ind(1)      
%       ind(2) = 1 + round((PopSize-Nelit-1)*rand) ;
%     end   
%     par(1,:) = NewPob(ind(1),:) ;
%     par(2,:) = NewPob(ind(2),:) ;
% IndCruz1=[];%zeros(2,numvars);
% IndCruz2=[];%zeros(2,numvars);
% 
% % if mod(igen,2)==0 %generacion par
% %   DesStan = std(NewPob) ;
% % j=1;
% %for ivar=1:numvars/2
% % for ivar=1:numvars/4
% % %  IndCruz1= [IndCruz1 par(1,j) par(2,j+1) ] ;
% % %  IndCruz2= [IndCruz2 par(2,j) par(1,j+1) ] ;
% %   IndCruz1= [IndCruz1 par(1,j) par(1,j+1) par(2,j+2) par(2,j+3) ] ;
% %   IndCruz2= [IndCruz2 par(2,j) par(2,j+1) par(1,j+2) par(1,j+3)] ;
% % %  IndCruz1= [IndCruz1 par(1,j) 0 par(2,j+2) 0 ] ;
% % %  IndCruz2= [IndCruz2 par(2,j) 0 par(1,j+2) 0 ] ;
% % %  PobCruz(2,:) = [PobCruz(2,:) par(2,j) par(1,j+1)] ;
% % %  [IndCruz1 ; IndCruz2]' ;
% %   j=j+4 ;
% % %  j=j+4 ;  
% % end�
% %    IndCruz1 = [mean(par(:,1))+DesStan(1)  mean(par(:,2))+DesStan(2)] ;
% %    IndCruz2 = [mean(par(:,1))-DesStan(1)  mean(par(:,2))-DesStan(2)] ;
% % else
% %    point = round(15.*sort(rand(2,1)))+ 1                       ;
% %    parcross = par(1,point(1,1):point(2,1))                     ;
% %    par(1,point(1,1):point(2,1)) = par(2,point(1,1):point(2,1)) ;
% %    par(2,point(1,1):point(2,1)) = parcross                     ;
% %    IndCruz1 = par(1,:);
% %    IndCruz2 = par(2,:);
% %    IndCruz1 = [par(1,1) par(2,2)] ;
% %    IndCruz2 = [par(2,1) par(1,2)] ;   
% % end
% 
%   esp = abs(par(:,2)-par(:,1)) ;
%   if mod(igen,2)==0 %generacion par
%     IndCruz1 = [par(1,2)-esp(1,1)  par(1,2)] ;
%     IndCruz2 = [par(2,2)-esp(2,1)  par(2,2)] ;
%     if (par(1,2)-esp(1,1)) < 0.0
%       IndCruz1 = [0.01  par(1,2)] ;
%     end
%     if (par(2,2)-esp(2,1)) < 0.0
%       IndCruz2 = [0.01  par(2,2)] ;       
%     end
%   else
%     IndCruz1 = [par(1,1)  par(1,1)+esp(2,1)] ;
%     IndCruz2 = [par(2,1)  par(2,1)+esp(1,1)] ;   
%   end  
% NewPob(ind(1),:)=IndCruz1 ;
% NewPob(ind(2),:)=IndCruz2 ;
% % PlotMesh(NewPob(ind(1),:))
% % PlotMesh(NewPob(ind(2),:)) 
%                         
% % PobCorr=NewPob(ind,:);
% %    [NewPobCorr] = CrossCorrectGAY(PobCorr,PopSize,numvars) ;
% % NewPob(ind,:)=NewPobCorr;
% 
% % PlotMesh(NewPob(ind(1),:))
% % PlotMesh(NewPob(ind(2),:))
%   end
% 
% 
% %%% MUTACION (Step mutation) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   for ind=1:PopSize-Nelit
%       if igen<=300
%         mut = (0.1 + 0.9*sin(pi/2*igen/300))*probmut ;
%       else
%         mut = probmut ;
%       end
%     if rand <= mut % se hace la mutaci�n con un 30% de probabilidad
%         rp = mean(NewPob(ind,:)) ;
%         esp_min = 1 ; % espesor minimo 1mm
%         ri_min  = 1 ; %radio interno minimo 1mm
%         esp = (esp_min + rand*((rp-ri_min)*2-esp_min));
% %         esp = (NewPob(ind,2)-NewPob(ind,1))*(0.1 + rand*1.9) ; %Aumenta o disminuye el espeso entre 0.1 y 2 veces
%         NewPob(ind,:) = [rp-esp/2  rp+esp/2] ;
% % %       coord  = (1 + round((numvars/2-1)*rand))*2-1  ;  %del 1 al 17
% %       coord  = (1 + round((numvars-1)*rand)) ;
% % %      change = NewPob(ind,coord) ;
% % %      while change == NewPob(ind,coord)
% % %        change = SecNorm(1 + round((size(SecNorm,1)-1)*rand)) ;
% % %      end
% %       jump= ;
% %       NewPob(ind,coord) = NewPob(ind,coord)*jump      ; 
% % %       NewPob(ind,coord+1) = NewPob(ind,coord+1)*jump ;
% % %      dir = round(rand) + 1 ; % dir = 1 limite inferior ; dir = 2 limite superior
% % %      limit = Bound(coord,dir)  ;
% % %      jump = rand ;
% % %      jump = round(1+size(SecNorm,1)*rand(PopSize,numvars))
% % %      NewPob(ind,coord) = jump*NewPob(ind,coord) + (1-jump)*limit ;
% 
% % PobCorr=NewPob(ind,:);
% %    [NewPobCorr] = CrossCorrectGAY(PobCorr,PopSize,numvars) ;
% % NewPob(ind,:)=NewPobCorr;
% 
%     end 
%   end    
% 
% %%% Nueva generaciones
%   NewPob = [NewPob ; Elit] ;
%   
%   for isol=1:PopSize
%     Sol(isol) = feval(funct,NewPob(isol,:),weight) ; 
%   end
%   graff=0;
%   if MinSol~=min(Sol)
%       graff=1;
%   end
%   MinSol = min(Sol);
%   PromSol = mean(Sol);
%   pos=find(Sol==MinSol);
% 
% %pause
% %figure(2)
% %%% Nuevas generaciones
% %plot3(NewPob(:,1),NewPob(:,2),Sol(:),'.')
% %contour(x,y,z);
% %grid
% %hold on
% %plot3(NewPob(pos(1),1),NewPob(pos(1),2),MinSol,'og')
% %plot(NewPob(:,1),NewPob(:,2),'.r')
% %hold off
% 
% igen2=igen2+1 ;
% BestInd = [BestInd ; NewPob(pos(1),:) MinSol] ;
% %MinSol
% fprintf(1,'---Generacion : %6i \n',igen) ;
% fprintf(1,'Mejor solucion: %6.8f \n',MinSol) ;
% ngen = [ngen igen2] ;
% minS = [ minS MinSol] ;
%    
% if lastMinSol~=MinSol
% %%% Graficos
% close(1)
% figure(1)
% %plot(igen,MinSol,'.g') %,igen,PromSol,'.b')
% lastMinSol=MinSol ;
% pos=find(BestInd(:,numvars+1)==MinSol) ;
% BestArea=BestInd(pos(1),1:1:numvars);
% subplot(2,1,1); plot(ngen,minS,'.g')
% xlim([-1 numgen*numit+1])
% xlabel('Generaciones')
% ylabel('Funcion')
% grid
% %close(1)
% % subplot(2,1,2); MEF2DBarraAG(BestArea)
% xlim([0 4])
% grid
% pause(0.1)
% end
% 
% % if (zcont == 0 & MinSol <= 0.187 )
% %    zcont = igen2 ;
% %    %pause
% % end
% end %bucle de igen
% 
% npass = round(PopSize*porpas/100) ;
% MaxSol = max(Sol) ;
% Sol    = Sol-MaxSol-epsilon  ;
% [Sol,ix] = sort(Sol,'descend') ; %'ascend'
% NewPobPass = [];
%  ieli = PopSize;
%   for inel = 1:npass
%     NewPobPass = [NewPobPass ; NewPob(ix(ieli),:)] ;
%     ieli = ieli-1;
%   end
% %   PopSize = PopSize-npass;
%   [InitPob] = RestPobInit((PopSize-npass),numvars);
%   NewPob = [NewPobPass ; InitPob] ;
% 
% 
% for isol=1:PopSize
%   Sol(isol) = feval(funct,NewPob(isol,:),weight) ; 
% end
%   
% % for indi=1:nfunc
% %   for isol=1:PopSize
% %     Func(indi,isol) = feval(funct,NewPob(isol,:),indi) ; %FuncExp(InitPob(isol,:)) ; %FuncTPN01a(InitPob(isol,:)) ;
% %   end
% %   media(indi) = mean(Func(indi,:));
% % end
% % weight=[0.1 1];
% % if nfunc==1
% %   Sol = Func(1,:) ;
% % elseif nfunc==2
% %   Sol = ((Func(1,:)./weight(1)).^2 + (Func(2,:)./weight(2)).^2).^0.5 ;
% % end
% 
% 
% 
% % figure(4)
% % for bind=1:size(NewPobPass,1)
% %     for indi=1:nfunc
% %         BestFunc(bind,indi) = feval(funct,NewPobPass(bind,1:1:numvars),[indi 0]) ; %FuncExp(InitPob(isol,:)) ; %FuncTPN01a(InitPob(isol,:)) ;
% %     end    
% % end
% % for indi=1:nfunc
% %    BestFunc2(1,indi) = feval(funct,BestInd(igen2,1:1:numvars),[indi 0]) ; %FuncExp(InitPob(isol,:)) ; %FuncTPN01a(InitPob(isol,:)) ;
% % end
% % plot(BestFunc(:,1),BestFunc(:,2),'.b',BestFunc2(1,1),BestFunc(1,2),'or')
% % grid
%   
% 
% hold off
% 
% 
% %%% Mejor soluci�n
% display('--- Best solution');
% MinSol = min(BestInd(:,numvars+1))
% display('--- Invidual');
% pos=find(BestInd(:,numvars+1)==MinSol) ;
% BestInd(pos(1),:)'
% feval(funct,BestInd(pos(1),1:1:numvars),[1 0])
% feval(funct,BestInd(pos(1),1:1:numvars),[0 1])
% 
% %plot3(BestInd(pos(1),1),BestInd(pos(1),2),BestInd(pos(1),3),'ob')
% %hold off
% % BestArea=BestInd(pos(1),1:1:numvars);
% % peso = sum(BestInd(pos(1),1:2:numvars)) ;
% % fprintf(1,'--- Valor proporcional al Peso Total: %6.8f \n',peso) ;
% % MEF2DBarraAG(BestArea)
% 
% % if nfunc==2
% % figure(3)
% % for bind=1:size(BestInd,1)
% %     for indi=1:nfunc
% %         BestFunc(bind,indi) = feval(funct,BestInd(bind,1:1:numvars),[indi 0]) ; %FuncExp(InitPob(isol,:)) ; %FuncTPN01a(InitPob(isol,:)) ;
% %     end    
% % end
% % plot(BestFunc(:,1),BestFunc(:,2),'o-')
% % grid
% % end