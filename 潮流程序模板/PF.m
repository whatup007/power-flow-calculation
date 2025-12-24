%*************************************
%功能：计算PF的主程序
%输入参数：FileName 算例文件名
%作者：
%模版作者：李佩杰
%模版生成时间：2023-11-12
%编写时间：
%调用方法：PF  case4gs
%        PF  case118
%*************************************
function converged = PF(FileName)
%clear all;

t0 = tic;
mpc = feval(FileName);
baseMVA = mpc.baseMVA;%基准容量
bus= mpc.bus;%节点信息
gen = mpc.gen;%发电机节点信息
branch = mpc.branch;%支路信息
%形成节点导纳矩阵
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

nb = size(bus, 1);
ng = size(gen, 1);
Cg = sparse(gen(:, 1), (1:ng)', gen(:, 8) > 0, nb, ng);  %% gen connection matrix
                                        %% element i, j is 1 if, generator j at bus i is ON
bus_gen_status = Cg * ones(ng, 1);      %% number of generators at each bus that are ON
on = find(gen(:, 8) > 0);      %% which generators are on?

%% 形成索引
ref = find(bus(:, 2) == 3 & bus_gen_status);   %% 平衡节点索引
pv  = find(bus(:, 2) == 2  & bus_gen_status);   %% PV节点索引
pq  = find(bus(:, 2) == 1 | ~bus_gen_status);   %% PQ 节点的索引
 %% initial state
converged = 0;
i = 0;
V0    = ones(size(bus, 1), 1);            %% 平启动
V = V0;%复数向量
Va = angle(V);%电压幅值
Vm = abs(V);  % 电压相角   


%% set up indexing for updating V
npv = length(pv);
npq = length(pq);


Sbus = Cg * (gen(on, 2) + 1j * gen(on, 3)) / baseMVA -(bus(:, 3)   + 1j * bus(:, 4) ) / baseMVA;

mis =               ;%计算不平衡分量，请补充
F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

%% 检查收敛精度
normF = norm(F, inf);
tol = 0.00001;%收敛精度
max_it = 50;%最大收敛次数

fprintf('\n it    max P & Q mismatch (p.u.)');
fprintf('\n----  ---------------------------');
fprintf('\n%3d        %10.3e', i, normF);
if normF < tol
    converged = 1;
    fprintf('\nConverged!\n');
end

%% 执行牛顿迭代
while (~converged && i < max_it)
    %%  更新迭代次数
    i = i + 1;
  
    %% 计算雅可比矩阵
   J= makeJacob(V,pv,pq,ref,Ybus);

    %% 计算修正量
    dx = J\(-F);


    %% 更新电压，请在下方补充 



        

    %% evalute F(x)
    mis =            ;%计算不平衡分量，请补充
    F = [   real(mis([pv; pq]));
        imag(mis(pq))   ];

    %% check for convergence
    normF = norm(F, inf);
    runtime = toc(t0);
    fprintf('\n%3d        %10.3e', i, normF);
    if normF < tol
        converged = 1;      
        fprintf('\nNewton''s method power flow (power balance, polar) converged in %d iterations.\n', i);
        printpf(runtime, converged, V,bus,branch,gen);%输出潮流结果
    end
end

if ~converged
   fprintf('\nNewton''s method power flow (power balance, polar) did not converge in %d iterations.\n', i);
end

