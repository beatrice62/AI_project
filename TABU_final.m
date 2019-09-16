%% TABU

%% 64-QAM, 2X2 @ SNR�� 39�� �������� �����Ǿ� ����
clc;  clear all;

%%����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 16; % modulation order (4�̻� �����ϵ��� ����, 4, 16, 64... ���)
nt= 2; %���׳� ��
iteration = 100000; % 10����
tabu_search_num = 30;   %% ������ M
num =1;
SNR_matrix = linspace(27,27,num);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���⵵ ����!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=log2(M); %% �� �ɺ��� ��Ʈ ��
number_of_bit = k* nt; %% �ѹ��� ������ ��Ʈ



hmod = modem.qammod('M', M, 'PHASEOFFSET', 0, 'SYMBOLORDER', 'Gray','INPUTTYPE', 'INTEGER');
hdemod = modem.qamdemod(hmod);

pt = modulate(hmod,[0:M-1]); % normalization ����
temp_Es = sum(real(pt).^2 + imag(pt).^2)/M*nt;
d = 1/sqrt(temp_Es);
pt = d*modulate(hmod,[0:M-1]); %normalization ����

for i=1:num
    
    error1(1,i) =0;
    
    for itn=1:iteration
        [itn];
        bit = randi([0,1],number_of_bit,1); % ��Ʈ ����
        symbol = reshape(bit,nt,number_of_bit/nt); % �ɺ� ����
        c= bi2de(symbol,'left-msb'); %msb�� ����
        x = pt(c+1); % mapping �Ϸ�
        x = transpose(x);
        
        SNR_dB = SNR_matrix(i);        % dB ����
        SNR_linear = 10^(SNR_dB/10);   % linear scale ����
        Es=1;                          % Es=1�� normalization �Ǿ� ����
        N0 = Es/SNR_linear;            % (N0 = Es/SNR [linear scale] ���� �̿�)
        sigma = (N0)^0.5;            % sigma^2 = N0/2 ���� �̿�
        noise = sigma /sqrt(2)* ( randn(nt, 1) + 1j * randn(nt, 1) );
        H= sqrt(1/2) * ( randn(nt, nt) + 1j * randn(nt, nt) );
        
        y_original = H*x+noise;
        
        
        
        
        %%%%%%%%%%%%%%%%%%% MMSE�� �ʱⰪ���� ����.%%%%%%%%%%%%%%%%
        Her_H = conj(transpose(H)); % H�� �츣�̽þ�
        G_MMSE = Her_H*( H*Her_H + (sigma)^2 * eye(nt) )^-1;
        % G_MMSE = inv(H);
        z1 = G_MMSE * y_original;
        
        result1=demodulate(hdemod,z1/d);
        initial_mmse = pt(result1+1);
        
        %%% real-valued system���� �ٲ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        H_real = [real(H), -imag(H) ; imag(H), real(H)];
        initial_sol = [real(initial_mmse).';imag(initial_mmse).']; %% Ÿ�� ��ġ���� ����� �ʱ���
        n_real = [real(noise);imag(noise)];
        x_real = [real(x);imag(x)];
        y_real = H_real*x_real + n_real;
        
        
        best_val = norm( y_real - H_real * initial_sol );
        best_sol = initial_sol;
        
        TS = [];
        TS = [TS initial_sol];
        tabu_number =[5];
        % T = 5*ones(2*nt*sqrt(M),2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %MLset=de2bi(0:M^nt-1,nt,M);
        %  ML_set_symbol= pt(MLset+1);
        %tabu_list_number = 5*ones(1,M^nt);
        %TS = zeros(2*nt,tabu_num);
        
        for index =1 : tabu_search_num
            
            
            %            TS_number = [tabu_num];
            
            candidate = neighbors(initial_sol,nt,M,d);
            %candidate = transpose(candidate);
            
            for integer =1 : 4*nt
                value(integer) = norm( y_real - H_real *  candidate(:,integer) );
            end
            
            [temp_mv,temp_mp]  = min(value);
            
            
            if temp_mv < best_val
                
                best_val = temp_mv;
                best_sol = candidate(:,temp_mp);
                initial_sol = best_sol;
            else
                
                for index_nn = 1: 4*nt
                    criteria=kron(candidate(:,index_nn),ones(1,size(TS,2)));
                    
                    is_there = (sum(abs(TS-criteria))==0);
                    temp_is_there = sum(is_there);
                    where = find(is_there==1);
                    
                    if temp_is_there >0 && tabu_number(where) > 0
                        value(index_nn) =1000;
                        
                    end
                end
                [temp_mv,temp_mp]  = min(value);
                initial_sol=candidate(:,temp_mp);
                
                criteria=kron(initial_sol,ones(1,size(TS,2)));
                is_there2 = (sum(abs(TS-criteria))==0);
                temp_is_there2 = sum(is_there2);
                where2 = find(is_there2==1);
                if temp_is_there2 >0
                    tabu_number(where2) = 5+1;
                else
                    
                    TS = [TS initial_sol];
                    tabu_number =[tabu_number 5+1];
                end
            end
            
            
            tabu_number = tabu_number - ones(1,size(tabu_number,2));
            
            
            
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        initial_sol = best_sol;
        initial_sol = reshape(initial_sol,nt,2);
        initial_sol = initial_sol(:,1) + 1j*initial_sol(:,2);
        
        result1=demodulate(hdemod,initial_sol/d);
        detection_bit1= de2bi(result1,k,'left-msb'); % msb�� ����
        detection_bit1 = reshape(detection_bit1, number_of_bit,1);
        
        difference1 = abs(bit - detection_bit1);      % ���� ��Ʈ�� ���� ��Ʈ�� ���Ͽ� ������ ã�´�
        error1(1,i) = error1(1,i)+ sum (difference1); % �� interation������ ������ ��������
        
        
        
        %%%% ������Ҵ��� ǥ��
        if mod(itn,10) == 0
            R = sprintf('present = %d, finally = %d', itn,iteration);
            disp(R)
        end
        
        
        
        
    end
    
    ber1 (1,i) = error1(1,i)/(number_of_bit*iteration);
    
end

figure(2)
SNR = 0:3:30;
semilogy(SNR_matrix,ber1,'b-+','linewidth',2);  hold on;



%2X2, 16QAM
ZF =[0.37175375	0.320545	0.260705	0.19783875	0.13834625	0.087967875	0.051948	0.028514	0.015006	0.00766	0.00392];
MMSE=[0.32816	0.27912625	0.2269075	0.17347375	0.1221725	0.078338875	0.04634725	0.025537625	0.013448625	0.006859	0.003538];
ML=[0.33787375	0.28646875	0.22963875	0.16845875	0.1088775	0.05836425	0.025494375	9.18E-03	2.83E-03	7.98E-04	2.01E-04];

semilogy(SNR,ZF,'r--','linewidth',1); hold on;
semilogy(SNR,MMSE,'b--','linewidth',1);
semilogy(SNR,ML,'g--','linewidth',1);

%
% xlim([0,26]);
% ylim([10^-4,1]);


grid on;




