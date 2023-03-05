function batch=lib_evaluation(batch, batch_cnt, signal_characteristics, param_method)

jumps=signal_characteristics.jumps;
N=signal_characteristics.N;
degree=signal_characteristics.degree;
y=signal_characteristics.y;

comparsion_breakpoints = find(diff(signal_characteristics.coeff(:,1))'~=0);
comparsion_breakpoints_1 = [comparsion_breakpoints, comparsion_breakpoints-1, comparsion_breakpoints+1];
comparsion_breakpoints_2 = [comparsion_breakpoints_1, comparsion_breakpoints-2, comparsion_breakpoints+2];

breakpoints=zeros(1,jumps);

for i=batch_cnt
    
    differences=sort(abs(batch{i}.diff_med_filt));
    differences_nonsorted=abs(batch{i}.diff_med_filt);
    position_differences=[];
    k=0;
    for j=1:jumps
        highest_value=differences(end-j+1);
        position_highest_value=find(differences_nonsorted==highest_value);
        if length(position_highest_value)>1
            len=length(position_highest_value);
            k=k+1;
            position_highest_value=position_highest_value(k);
            if k==len
                k=0;
            end
            
        end
        
        if position_highest_value>2 && position_highest_value<(N-2)
        differences_nonsorted(position_highest_value+1)=0;
        differences_nonsorted(position_highest_value+2)=0;
        differences_nonsorted(position_highest_value-1)=0;
        differences_nonsorted(position_highest_value-2)=0;
        
            if find((position_highest_value+1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value+1];
            end
            if find((position_highest_value+2)~=position_differences)==0
                position_differences=[position_differences,position_highest_value+2];
            end
            if find((position_highest_value-1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value-1];
            end
            if find((position_highest_value-2)~=position_differences)==0
                position_differences=[position_differences,position_highest_value-2];
            end
        end
        
        if position_highest_value==(N-2)
        differences_nonsorted(position_highest_value+1)=0;
        differences_nonsorted(position_highest_value-1)=0;
        differences_nonsorted(position_highest_value-2)=0;
        
            if find((position_highest_value+1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value+1];
            end
            if find((position_highest_value-1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value-1];
            end
            if find((position_highest_value-2)~=position_differences)==0
                position_differences=[position_differences,position_highest_value-2];
            end
        end
        
        if position_highest_value==(N-1)
        differences_nonsorted(position_highest_value-1)=0;
        differences_nonsorted(position_highest_value-2)=0;
        
           
            if find((position_highest_value-1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value-1];
            end
            if find((position_highest_value-2)~=position_differences)==0
                position_differences=[position_differences,position_highest_value-2];
            end
        end
        
        if position_highest_value==2
        differences_nonsorted(position_highest_value-1)=0;
        differences_nonsorted(position_highest_value+1)=0;
        differences_nonsorted(position_highest_value+2)=0;
        
            if find((position_highest_value+1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value+1];
            end
            if find((position_highest_value+2)~=position_differences)==0
                position_differences=[position_differences,position_highest_value+2];
            end
            if find((position_highest_value-1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value-1];
            end
           
        end
        
        if position_highest_value==1
        differences_nonsorted(position_highest_value+1)=0;
        differences_nonsorted(position_highest_value+2)=0;
        
            if find((position_highest_value+1)~=position_differences)==0
                position_differences=[position_differences,position_highest_value+1];
            end
            if find((position_highest_value+2)~=position_differences)==0
                position_differences=[position_differences,position_highest_value+2];
            end
            
        end
        
        differences=sort(differences_nonsorted);
        highest(j)=highest_value;
        breakpoints(j)=position_highest_value;
        
    end
    
    pom=length(position_differences);

    rest=differences(1+pom:end-5);
    average_rest=mean(rest);
    average_highest=mean(highest);
   

    NoB=0;
    for j=1:jumps
%         breakpoints(j)=find(differences_nonsorted==highest(j)); %TODO: upravit
        if sum(comparsion_breakpoints_2==breakpoints(j))>0
            NoB=NoB+1;
        end
    end
    
    AAR=average_highest/average_rest;
    MMR=min(highest)/max(rest);
    
    disp(['Statistics for the ', batch{batch_cnt}.algorithm, ' algortihm:'])
    disp([num2str(jumps), ' highest breakpoints in reconstructed signal: ', num2str(sort(breakpoints))])
    disp(['Statistics: AAR: ', num2str(AAR), ', MMR: ', num2str(MMR), ', NoB: ', num2str(NoB) ])
    
    % MSE reconstructed signal from used reconstruction method
    MSE_recon_from_method = (1/N)*norm(signal_characteristics.signal_clean-batch{batch_cnt}.y_recon)^2;
    y_error=batch{batch_cnt}.y_recon-signal_characteristics.signal_clean;
    SNR_dB_method = 10*log10(sum(signal_characteristics.signal_clean.^2)/sum(y_error.^2));
    % MSE reconstructed signal from 5 highest values
    % refitting of each detected segment by using least squares method
    recon_highest=zeros(N,degree);
    breakpoints_5_highest_sort=[0,sort(breakpoints),N];
    for ii=1:length(breakpoints_5_highest_sort)-1
        start=breakpoints_5_highest_sort(ii)+1; %start position of segment
        stop=breakpoints_5_highest_sort(ii+1); %stop position of segment

        % getting new paramtererization coefficients
        beta_det=lib_least_squares_fit(start, stop, y(start:stop), param_method);

        % creating new parameterizaion coefficients
        for jj=1:degree
            recon_highest(start:stop,jj)=ones(stop-start+1,1).*beta_det(jj);
        end
    end

    signal_recon_breakpoints_5_highest = lib_recon_from_params(reshape(recon_highest,degree*N,1),N,degree, param_method); % reconstructed signal
    MSE_denoised_recon_5_highest = (1/N)*norm(signal_characteristics.signal_clean-signal_recon_breakpoints_5_highest)^2;
    y_error=signal_recon_breakpoints_5_highest-signal_characteristics.signal_clean;
    SNR_dB_recon_highest = 10*log10(sum(signal_characteristics.signal_clean.^2)/sum(y_error.^2));
    
    batch{batch_cnt}.signal_recon_highest=signal_recon_breakpoints_5_highest;
    
    % Display signal
    figure ('Name', 'Reconstructed signal of highest breakpoins', 'NumberTitle', 'off'); 
    title ('Reconstructed signal of highest breakpoins')
    plot(signal_characteristics.signal_clean, 'y');
    hold on
    plot(signal_characteristics.y, 'b');
    
    plot(batch{batch_cnt}.y_recon, 'r');
    
    plot(batch{batch_cnt}.signal_recon, 'g');
    
    plot(signal_recon_breakpoints_5_highest, 'c');
    legend('clean','noisy',['recon ', batch{batch_cnt}.algorithm], 'recon from detected BP', 'recon from highest BP');
    
    
    MSE_breakpoints_recon = (1/N)*norm(signal_characteristics.signal_clean-batch{batch_cnt}.signal_recon)^2;
    y_error=batch{batch_cnt}.signal_recon-signal_characteristics.signal_clean;
    SNR_dB_recon_from_param = 10*log10(sum(signal_characteristics.signal_clean.^2)/sum(y_error.^2));
    
    
    disp(['Statistics: MSE of signal reconstructed from ', batch{batch_cnt}.algorithm ,': ', num2str(MSE_recon_from_method)])
    disp(['Statistics: MSE of signal reconstructed from detected breakpoint ', batch{batch_cnt}.algorithm ,': ', num2str(MSE_breakpoints_recon)])
    disp(['Statistics: MSE of signal reconstructed from ', num2str(jumps) ,' highest breakpoints: ', num2str(MSE_denoised_recon_5_highest)])
    disp(['Statistics: SNR of signal reconstructed from ', batch{batch_cnt}.algorithm ,': ', num2str(SNR_dB_method)])
    disp(['Statistics: SNR of signal reconstructed from detected breakpoint ', batch{batch_cnt}.algorithm ,': ', num2str(SNR_dB_recon_from_param)])
    disp(['Statistics: SNR of signal reconstructed from ', num2str(jumps) ,' highest breakpoints: ', num2str(SNR_dB_recon_highest)])
    
    %     
    legend_evaluation=[{'AAR'},{'MMR'},{'NoB'},{'MSE_of_signal_recon_from_method'}...
        ,{'MSE_of_signal_recon_from_detected_param'},{'MSE_of_signal_recon_from_highest_breakpoints'}...
        ,{'SNR_of_signal_recon_from_method'},{'SNR_of_signal_recon_from_detected_param'}...
        ,{'SNR_of_signal_recon_from_highest_breakpoints'}];
    batch{batch_cnt}.legend_evaluation=legend_evaluation;
    
    results_of_evaluation=[AAR,MMR,NoB,MSE_recon_from_method,...
        MSE_breakpoints_recon,MSE_denoised_recon_5_highest,...
        SNR_dB_method,SNR_dB_recon_from_param...
        SNR_dB_recon_highest];
    
    batch{batch_cnt}.evaluation_results=results_of_evaluation;
end