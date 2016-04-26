clear all; close all, clc

%% PARAMETRY OBIEKTU i ESTYMACJI
iter= 4000;      % liczba iteracji
t= [1:iter];     % wektor iteracji
dwe= 1;          % odchylenie standardowe wejúcia (szum normalny)
dsz= 0.01;    % odchylenie standardowe szumu pomiarowego

%Wybierz algorytm: 1-zwyk≥y WRLS, 2-WRLS(zmienna lambda) pierwszy 3-WRLS(zmienna lambda) drugi
algorytm = 3;

switch(algorytm)


    case(3)
        
        lambda =1;       % warunek poczπtkowy lambdy
        % wspÛ≥czynniki czu≥oúci 
        L=2;
        Ka=6;
        Kb=12;
        alpha=1-1/(Ka*L)
        beta=1-1/(Kb*L)
        gamma=1.01;
        psi=10e-8;
        
        vare=0;
        varv=0;
        qkr=10;
        
	otherwise
    	return;
end

% Parametry obiektu =a(t)*x +b(t)
x1 = ones(1,300); x2 = zeros(1,300); x3 = linspace(0,2,400); x4 = linspace(2,-2,200);
x5 = linspace(-2,0,200); x6 = zeros(1,300); x7 = ones(1,300);
a=[x1 x2 x3 x4 x5 x6 x7 x1 x2 x3 x4 x5 x6 x7];
% a=[ones(1,500) zeros(1,500) ones(1,500) zeros(1,500)]; 

% b=[zeros(1,1000) ones(1,1000)];
b=[ones(1,1000) zeros(1,1000) ones(1,1000) zeros(1,1000)];

A=[];           % macierz pod estymaty parametru a - historia
B=[];           % macierz pod estymaty parametru b - historia
LAMBDA=[];      % macierz pod parametr lambda - historia
SIG_OUT=[];     % macierz pod b≥πd predykcji

%% generacja sygna≥u testowego i odpowiedzi
sig_test=dwe*randn(1,iter);                         %sygna≥ testowy (szum o odchyleniu dwe)
sig_out_ideal=sig_test.*a+b;                        %idealna odpowiedü
sig_out_real=sig_test.*a+b+dsz*randn(1,iter);       %pomiar (idealna odpowiedü + szum pomiarowy(odchylenie dsz))


%% estymacja WRLS 

%warunki poczatkowe
est=[0;0];              %wg. instrukcji na Êwiczenia z RLS
P=eye(2).*10e3;

for i=1:iter
    
    switch(algorytm)
	case(1)
            
	case(2)
        
        L= -(2*(sig_out_real(i)-[sig_test(i); 1]'*est).^2); % modyfikacja 
        lambda= lambdamin+(1-lambdamin).*(2.^L);                   
        lambda= min(lambda,1);

    case(3)
        
        q=[sig_test(i);1]'*P*[sig_test(i);1];
        qkr=alpha*qkr+(1-alpha)*q;
        
        e=sig_out_real(i)-[sig_test(i); 1]'*est;
    
        vare=alpha*vare+(1-alpha)*e^2;
        varv=beta*varv+(1-beta)*e^2;
    
        if sqrt(vare)<=(gamma*sqrt(varv));  
            lambda=1
        else
            lambda=min(qkr*varv/(psi+abs(vare-varv)),1);
        end
    otherwise
    	return;
    end
   
    c=inv([sig_test(i); 1]'*P*[sig_test(i); 1]+lambda);
    K=P*[sig_test(i); 1]*c;                   % Korekta
                                                      
    P=1/lambda .*(P-K*[sig_test(i); 1]'*P);
    e=sig_out_real(i)-[sig_test(i); 1]'*est;  %b≥πd predykcji
    
    est=est+K*(e);
  
    
    A=[A est(1)];
    B=[B est(2)];
    LAMBDA = [LAMBDA lambda];
    SIG_OUT =  [SIG_OUT e];
  
end
% B≥Ídy estymacji
delta_a1=a-A;
delta_b1=b-B;

SIG_OUT = SIG_OUT.^2;

%% prezentacja wynikÛw
figure
subplot(311);
plot(t, A,t, a,'r'); grid on; title('Estymacja parametru A');...
legend('estymata parametru a','parametr a'); grid on;hold on; ylabel('Amplituda'); 

subplot(312);
plot(t, B,t, b,'r'); grid on; title('Estymacja parametru B');...
legend('estymata parametru b','parametr b'); grid on;hold on; ylabel('Amplituda'); 

subplot(313);

plot(t, LAMBDA); grid on;
legend('parametr lambda'); grid on;hold on; ylabel('lambda'); title('Zmiana parametru lambda');

% subplot(514);
% plot(t, delta_a1,'g', t, delta_b1, 'r'); grid on;
% title('B≥πd estymacji parametrÛw'); grid on;hold on; ylabel('Amplituda'); legend('delta a','delta b');
% 
% subplot(515);
% plot(t, SIG_OUT ,'g'); grid on;
% title('Kwadrat b≥Ídu predykcji odpowiedzi obiektu'); grid on; hold on; ylabel('Amplituda'); legend('Kwadrat b≥Ídu predykcji');
