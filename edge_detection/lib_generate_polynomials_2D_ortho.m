function [P_big] = lib_generate_polynomials_2D(signal_characteristics,param_generation_signal)
% This function generates 2D polynomials

if ~isfield(param_generation_signal, 'type'), param_generation_signal.type = 4 ; end 
if ~isfield(param_generation_signal, 'plot_polynomials'), param_generation_signal.plot_polynomials = 0; end 

N = signal_characteristics.N;
M = signal_characteristics.M;
degree = signal_characteristics.degree;

type = param_generation_signal.type;

% type of polynomials
% 1 - standard 
% 2 - Legendre
% 3 - Chebyshev
% 4 - random orthogonal
%% type of polynomials
switch type
    case 1
        % original
        p{1}=@(x) x.^0;
        p{2}=@(x) x.^1;
        p{3}=@(x) x.^2; 
        
        lim_min=0;
        lim_max=1;
    case 2
        % Legendere
        p{1} = @(x) x.^0;
        p{2} = @(x) x;
        p{3} = @(x) 3/2*x.^2 - 1/2;
        p{4} = @(x) 5/2*x.^3 - 3/2*x;
        p{5} = @(x) 35/8*x.^4 - 30/8*x.^2 + 3/8;
        p{6} = @(x) 63/8*x.^5 - 70/8*x.^3 + 15/8*x;
        p{7} = @(x) 231/16*x.^6 - 315/16*x.^4 + 105/16*x.^2 - 5/16;

        lim_min=0;
        lim_max=1;

    case 3
        %Chebyshev
        p{1} = @(x) x.^0;
        p{2} = @(x) x;
        p{3} = @(x) 2*x.^2 - 1;

        lim_min=-1;
        lim_max=1;
        
end


%% generating of polynomials of length N
switch type
    case {1,2,3}
        x_N=linspace(lim_min,lim_max,N);



        P_N=ones(N,degree);
        for i=1:degree
            P_N(:,i)=p{i}(x_N);
        end
        
    case 4 
        A=zeros(N,degree);
        timevec=(1:N)/N;

        for i=1:degree
            d=timevec.^(i-1);
            A(:,i)=d';
            norms(i) = norm(A(:,i));
            A_norm(:,i)=d'/norm(d);
            norm(A_norm(:,i)); %should be one
        end

        B=randn(degree, degree);

        U_new = A_norm * B * diag( 2*rand(3,1)-1);  %make linear combination and randomize lengths

        for i=1:degree
            U_new(:,i)= U_new(:,i)./norm(U_new(:,i));
            norm(U_new(:,i))
        end

        [P_N,~,~] = svd(U_new,'econ');
            
end

%% generating of polynomials of length M
switch type
    case {1,2,3}
        x_M=linspace(lim_min,lim_max,M);



        P_M=ones(M,degree);
        for i=1:degree
            P_M(:,i)=p{i}(x_M);
        end
        
    case 4
        
        A=zeros(M,degree);
        A_norm=A;
        timevec=(1:M)/M;

        for i=1:degree
            d=timevec.^(i-1);
            A(:,i)=d';
            norms(i) = norm(A(:,i));
            A_norm(:,i)=d'/norm(d);
            norm(A_norm(:,i))  %should be one
        end

        B=randn(degree, degree);

        U_new = A_norm * B * diag( 2*rand(3,1)-1);  %make linear combination and randomize lengths

        for i=1:degree
            U_new(:,i)= U_new(:,i)./norm(U_new(:,i));
            norm(U_new(:,i))
        end

        [P_M,~,~] = svd(U_new,'econ');
            
end

%% generating of 2D polynomials

P_big=zeros(N*M,degree*degree);
pom=0;
for n=1:degree
    for m=1:degree
        pom=pom+1;
        P_big(:,pom)=kron(P_N(:,n),P_M(:,m));
    end
end

%% plot 2D polynomials
% if param_generation_signal.plot_polynomials == 1
%     for i=1:pom 
%         figure ('Name', ['Obrazek', num2str(i)])
%         image (reshape(P_big(:,i),M,N), 'CDataMapping','scaled')
%         colormap (gray)
%         colorbar
% %         imwrite(uint16(reshape(P_big(:,i),M,N)),['fig_polynomial_2D_ortho_' num2str(i) '.png'])
%    
%     end
% end