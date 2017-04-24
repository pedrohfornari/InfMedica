% Universidade Federal de Santa Catarina
% Centro Tecnologico
% Departamento de Engenharia Eletrica e Eletronica
% Disciplina: Introducao a Informatica Medica (EEL7307)
% Professora: Christine Fredel Boos
% Semestre letivo: 2017/1

% Atividade pratica da aula de Redes Neurais Artificiais

close all; clear all; clc

local_rna = '~/Documents/InfMedica/03.Neural_Networks/rna_mlp/'; % pasta onde estao/ficarao os arquivos

% Padroes de entrada e saida para o treinamento da rede
% load([local_rna,'pesos_pdrs.mat']);
arq_dados_rna = [local_rna,'Padroes e pesos.xls'];  % endereco da planilha do Excel
x = xlsread(arq_dados_rna,'padroes','B22:F119');    % padroes de entrada
dk = xlsread(arq_dados_rna,'padroes','R22:V26');    % padroes de saida

% Pesos sinopticos iniciais para a camada oculta e de saida
w = xlsread(arq_dados_rna,'pesos'); % pesos iniciais
wji = w(5:103 , 2:15); % pesos iniciais (camada oculta)
wkj = w(5:19 , 18:22); % pesos iniciais (camada de saida)
qtde_padrao = size(x,2); % quantidade de padroes de treinamento

% Parametros de configuracao da rede e do treinamento
n1 = ;                                          % neuronios na camada de entrada
n2 = ;                                          % neuronios na camada oculta
n3 = ;                                          % neuronios na camada de saida
tx_lrn = ;                                      % taxa de aprendizagem
coef_m = ;                                      % coeficiente de momento
erro_min_epc = ;                                % erro minimo por epoca
max_epocas = ;                                  % maximo numero de epocas
bias_1 = ;                                      % bias da camada oculta
bias_2 = ;                                      % bias da camada de saida

% Inclinacao das funcoes de ativacao
k1 = ;                                          % inclinacao da funcao na camada oculta
k2 = ;                                          % inclinacao da funcao na camada de saida

% Inicializacao de variaveis que armazenarao os erros da rede e dos padroes
erro_pdr = zeros(max_epocas,qtde_padrao);
erro_epc = zeros(max_epocas,1);
emq_epc = zeros(max_epocas,1);

% Grafico para visualizacao do erro de treinamento
grafico_rna = figure('NumberTitle','off','Color',[.94 .94 .94],...
    'Name','Treinamento da Rede Neural Artificial - grafico do erro da rede',...
    'Units','normalized', 'Position',[0.2556 0.2591 0.4392 0.5208]);
plot(0,0,'tag','dados_rna');
set(gca, 'FontSize', 8 , 'Position', [0.09 0.2 0.885 0.75]);
xlabel('Epocas de treinamento','FontSize',8);
ylabel('Valor do erro','FontSize',8);
h_cancel = uicontrol('Style','pushbutton', 'String','CANCELAR',...
    'FontUnits','normalized', 'FontSize',0.35, 'FontName','Verdana',...
    'ForegroundColor',[0.847 0.161 0],'FontWeight','bold',...
    'Units','normalized', 'Position',[0.435 0.01 0.14 0.075],'Callback', 'cancelar = 1;');
align(h_cancel,'Center','none');
cancelar = 0;
movegui(grafico_rna,'center')

% Inicio do treinamento da rede
for epc = 1:max_epocas                          % p/ cada epoca de treinamento da rede
    for pdr = 1:qtde_padrao                     % p/ cada padrao apresentado para a rede
        
        xi = [bias_1*ones(1,qtde_padrao);x];    % valores de entrada da rede
        
        % Ativacao dos neuronios da camada intermediaria (oculta)
        vj_aux = zeros(size(xi,1),n2);
        for j = 1:n2;
            vj_aux(:,j) = wji(:,j).*xi(:,pdr);
        end
        vj = sum(vj_aux)';                      % ativacao dos neuronios ocultos
        
        % Saida (yj) dos neuronios da camada intermediaria (oculta) -> equacao 1
        
        
        % Ativacao (vk) dos neuronios da camada de saida -> equacao 2
        
        
        % Saida (yk) dos neuronios da camada de saida -> equacao 3
        
        
        % Calculo dos erros local (ek) e do padrao (erro_pdr)
        % erro local -> erro de cada neuronio de saida
        ek = dk(:,pdr) - yk;
        % erro medio quadratico por padrao (considerando todos os neuronios de saida)
        erro_pdr(epc,pdr) = sum((ek.^2)/2);                                     
                
        % Calculo dos erros da camada de saida (delta_k) -> equacao 4
        
        
        % calculo dos erros da camada oculta (delta_j), considerando que a funcao de
        % ativacao da a funcao logistica
        delta_j = (yj.*(1-yj)) .* (wkj(2:end,:) * delta_k);
        
        % Acumulo do ajuste dos pesos sinopticos (dwkj e dwji)
        if pdr == 1;
           dwkj = tx_lrn .* (yj_aux * delta_k');
           dwji = tx_lrn .* (xi(:,pdr) * delta_j');
        else
           dwkj = dwkj + (tx_lrn .* (yj_aux * delta_k'));
           dwji = dwji + (tx_lrn .* (xi(:,pdr) * delta_j'));
        end
    end
    
    % C�lculo do erro da �poca (erro_epc)
    erro_epc(epc,1) = sum(erro_pdr(epc,:));      % erro da �poca
    emq_epc(epc,1) = mean(erro_pdr(epc,:));      % erro m�dio da �poca
    
    % Vizualiza��o do gr�fico com o erro de treinamento por �poca e por padr�o
    set(findobj('tag','dados_rna'),'XData', (1:1:epc),'YData',erro_epc(1:epc,1));
    drawnow; pause(.05)
    if cancelar == 1
        break;
    end
    % Verifica��o de crit�rio de parada: EMQ m�nimo por �poca
    % se o EMQ da �poca � menor que erro_min_epc -> FIM DO TREINAMENTO
    if ( emq_epc(epc,1) < erro_min_epc );
        break;
    end
    
    % C�lculo dos novos pesos sin�pticos (wkj_novo e wji_novo) -> equa��o 5
    if epc == 1
        
    else
        
    end
    
    % Atualiza��o dos pesos sin�pticos
    wkj_antigo = wkj;                           % pesos da camada oculta (t-1)
    wji_antigo = wji;                           % pesos da camada de sa�da (t-1)
    wkj = wkj_novo;                             % pesos da camada oculta (t)
    wji = wji_novo;                             % pesos da camada de sa�da (t)
end

vetor_epc = get(findobj('tag','dados_rna'),'XData');
vetor_erro = get(findobj('tag','dados_rna'),'YData');
close(grafico_rna);

figure('NumberTitle','off','Color',[.94 .94 .94],...
    'Name','Treinamento da RNA - gr�fico do erro da rede',...
    'Units','normalized', 'Position',[0.2556 0.2591 0.4392 0.5208]);
plot(vetor_epc, vetor_erro, '-b', 'LineWidth',1.5)
set(gca, 'FontSize', 8, 'Units','normalized', 'Position',[0.0715 0.11 0.9 0.8625]);
text(max(get(gca,'XLim')/2),max(get(gca,'YLim')/2),...
    ['Erro de treinamento = ',num2str(vetor_erro(end))],...
    'HorizontalAlignment','center','FontSize',9,'FontName','Verdana',...
    'Color',[0.847 0.161 0])
xlabel('�pocas de treinamento','FontSize',9);
ylabel('Valor do erro','FontSize',9);

% Final do treinamento
salvar = questdlg('Deseja salvar a rede?','RNA', 'Sim','N�o','Sim');
if strcmp(salvar,'Sim') == 1
    arq_rna = [local_rna,'Rede neural.xls'];    % planilha com o resultado do treinamento
    [stat1,msg1] = xlswrite(arq_rna,wji,'Teste da rede','T5:AG103');
    [stat2,msg2] = xlswrite(arq_rna,wkj,'Teste da rede','BE5:BI19');
    [stat3,msg3] = xlswrite(arq_rna,[vetor_epc', emq_epc(1:epc,:), ...
        erro_pdr(1:epc,:)], 'Aprendizagem da rede','A2');
    if ( stat1 == 0 ) || ( stat2 == 0 ) || ( stat3 == 0 )
        disp(' '); disp('Exportar pesos ocultos:'); disp(msg1);
        disp(' '); disp('Exportar pesos de sa�da:'); disp(msg2);
        disp(' '); disp('Exportar erros epc e pdr:'); disp(msg3); disp(' ');
        errordlg('Verifique o janela de comando','Erro ao exportar! Tente novamente...');
    else
        clc;
        msgbox({'   Rede pronta para teste!';''},'RNA');
    end
end
clearvars -except n1 n2 n3 bias_1 bias_2 phi1 phi2 k1 k2 wji wkj