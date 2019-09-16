%% GENETIC

%% ����ƽ �����ϼ�(5/4) ����Ʈ�� ��Ʈ�Ͻ�(���� ����), uniqueness �̿�
% (�ߺ� ��� O, �귿 �� ������ ��� 0)

clc; close all; clear all;
%% ��ⷹ�̼� order, ���׳� ����, iteration Ƚ�� ����!!!!!
M = 16; % modulation order (4�̻� �����ϵ��� ����, 4, 16, 64... ���)
nt= 2;  % ���׳� ��
iteration = 10000; 

%% �׷��� ����!!!!
num = 3;
starting = 15;
ending = 21;

%% genetic parameter ����!!!
population = 30;
generation = 15;
elite = 3;
selection_probability = 0.1;  % Ȯ�� ���� 1�ǵ���
crossover_probability = 0.8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mutation_probability = 1 - selection_probability - crossover_probability;



%%
SNR_matrix = linspace(starting,ending,num);

k=log2(M); %% �� �ɺ��� ��Ʈ ��
number_of_bit = k* nt; %% �ѹ��� ������ ��Ʈ

hmod = modem.qammod('M', M, 'PHASEOFFSET', 0, 'SYMBOLORDER', 'Gray','INPUTTYPE', 'INTEGER');
hdemod = modem.qamdemod(hmod);

pt = modulate(hmod,[0:M-1]); % normalization ����
temp_Es = sum(real(pt).^2 + imag(pt).^2)/M*nt;
d = 1/sqrt(temp_Es);
pt = d*modulate(hmod,[0:M-1]); %normalization ����
pt_real = -( sqrt(M)-1 ):2:  sqrt(M)-1 ;
pt_real = d* pt_real;

for i=1:num
    
    error1(1,i) =0;
    uniqueness(1,i)=0;
    for itn=1:iteration
        
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
        
        
        %best_val = norm( y_real - H_real * initial_sol );
        %best_sol = initial_sol;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% ����ƽ
        
        
        random_solution=pt_real(randi([1 sqrt(M)],2*nt,population-1));
        current_generation = [initial_sol random_solution]; %% �ʱ� generation ���� �Ϸ�
        for indexx = 1: population
            current_norm(indexx) =  norm( y_real - H_real * current_generation(:,indexx) );
            weight(indexx) = 1/current_norm(indexx);
        end
        
        
        all_generation = current_generation;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for number_generation = 1:generation-1
            
            %%% ����Ʈ�� ���� �����
            
            [~, position] = sort(current_norm);
            next_generation = current_generation(:,position(1:elite)); %% ����Ʈ �ֱ� �Ϸ�
            
            
            for indexss = 1:elite
                next_norm(indexss) = norm( y_real - H_real * next_generation(:,indexss) );
            end
            
            
            dice = rand(1, population-elite);
            temp_selection = sum(dice<=selection_probability);
            temp_crossover = sum(dice > selection_probability & dice < selection_probability + crossover_probability);
            temp_mutation = population - elite - temp_selection - temp_crossover;
            
            %%% ����Ʈ �̿ܿ� ���� �����
            for sel = 1: temp_selection
                
                yes_overlap=1;
                sel_weight = weight;
                while yes_overlap==1 % �ߺ��̸� ��� �ݺ�
                    random_number = fortune_wheel(sel_weight);
                    temp_vector = current_generation(:,random_number);
                    temp_norm = norm(y_real - H_real * temp_vector);
                    yes_overlap = overlap_confirm(temp_norm, next_norm);
                    sel_weight(random_number) = sel_weight(random_number)*0.5; % ��ġ�� ���� Ȯ�� ����
                end
                next_generation = [next_generation      temp_vector];
                next_norm = [next_norm  temp_norm];
            end
            
            
            for cross = 1: temp_crossover
                yes_overlap=1;
                cross_weight = weight;
                while yes_overlap==1
                    random_number1 = fortune_wheel(cross_weight);
                    random_number2 = fortune_wheel(cross_weight);
                    chromosome1 = current_generation(:,random_number1);
                    chromosome2 = current_generation(:,random_number2);
                    new_chromosome=[chromosome1(1:nt); chromosome2(nt+1:end)];
                    temp_norm = norm(y_real - H_real * new_chromosome);
                    yes_overlap = overlap_confirm(temp_norm, next_norm);
                    cross_weight(random_number1) = cross_weight(random_number1)*0.5;
                    cross_weight(random_number2) = cross_weight(random_number2)*0.5;
                end
                next_generation = [next_generation     new_chromosome];
                next_norm = [next_norm  temp_norm];
            end
            
            
            for mut = 1: temp_mutation
                yes_overlap=1;
                mut_weight = weight;
                while yes_overlap==1
                    random_number = fortune_wheel(mut_weight);
                    mutation_chromosome = current_generation(:,random_number); % �����ϰ� �ϳ� ����
                    mutation_chromosome(randperm(2*nt,1),1) = pt_real(randperm(sqrt(M),1)); % ���� �߻�
                    temp_norm = norm(y_real - H_real * mutation_chromosome);
                    yes_overlap = overlap_confirm(temp_norm, next_norm);
                    mut_weight(random_number) = mut_weight(random_number)*0.5;
                end
                next_generation = [next_generation     mutation_chromosome];
                next_norm = [next_norm  temp_norm];
            end
            
            
            current_generation = next_generation;
            current_norm = next_norm;
            next_norm = [];
            next_generation = [];
            
            for indexx = 1: population
                weight(indexx) = 1/current_norm(indexx);
            end
            
            all_generation = [all_generation current_generation];
            
            
        end
        
        [u_value,u_index]=unique(all_generation.','rows','stable');
        uniqueness(1,i) =  uniqueness(1,i) + length(u_index);
        
        
        
        [~, position] = sort(current_norm);
        best_sol = current_generation(:,position(1));
        
        
        
        
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
            R = sprintf('present = %d, finally = %d, snr=%d, finally=%d', itn,iteration,i,num);
            disp(R)
        end
        
        
        
        
    end
    
    ber1 (1,i) = error1(1,i)/(number_of_bit*iteration);
    average_uniqueness(1,i) = uniqueness(1,i)/iteration;
end


semilogy(SNR_matrix,ber1,'b-+','linewidth',2);  hold on;

% 
SNR =0:3:30;
%2X2
ZF =[0.37175375	0.320545	0.260705	0.19783875	0.13834625	0.087967875	0.051948	0.028514	0.015006	0.00766	0.00392];
MMSE=[0.32816	0.27912625	0.2269075	0.17347375	0.1221725	0.078338875	0.04634725	0.025537625	0.013448625	0.006859	0.003538];
ML=[0.33787375	0.28646875	0.22963875	0.16845875	0.1088775	0.05836425	0.025494375	9.18E-03	2.83E-03	7.98E-04	2.01E-04];

% 4X4
% ZF =[0.412235	0.372058125	0.320550625	0.26050375	0.19843625	0.138708813	0.088679813	0.051768438	0.028429438	0.014968188	0.007778188];
% MMSE=[0.33512	0.291685625	0.246561875	0.200385	0.15353875	0.108699063	0.070181438	0.04134475	0.022762125	0.012040625	0.006276438];
% ML=[0.34476875	0.296265625	0.24349625	0.188568125	0.125351875	0.05893225	0.015814563	2.26E-03	2.04E-04	1.61E-05	2.38E-06];


%%%64QAM 2x2 (������ ���)

% SNR = 0:2:40;
% ZF =[	0.4155425	0.393333333	0.3652825	0.3340225	0.299315833	0.261959167	0.224044167	0.186093333	0.149173333	0.117381667	0.086940833	0.06284	0.0435275	0.029424167	0.0195275	0.01284	0.007963333	0.005235	0.003265	0.002036667	0.001355];
% MMSE	=[0.387459167	0.3638625	0.336806667	0.30858	0.277144167	0.243766667	0.2098225	0.176164167	0.141248333	0.111436667	0.083216667	0.060094167	0.041701667	0.0283325	0.018766667	0.0123575	0.007616667	0.005055	0.003130833	0.001965	0.001310833];
% ML	=[0.395018333	0.372126667	0.344130833	0.313820833	0.281000833	0.245559167	0.20998	0.17409	0.137885	0.104744167	0.0724375	0.046720833	0.027826667	0.015083333	0.007610833	0.003715833	0.001670833	0.000774167	0.0003375	0.000129167	5.08E-05];
% 


semilogy(SNR,ZF,'r--','linewidth',1); hold on;
semilogy(SNR,MMSE,'b--','linewidth',1);
semilogy(SNR,ML,'g--','linewidth',1);

%
% xlim([0,26]);
% ylim([10^-4,1]);


grid on;




