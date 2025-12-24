"""
牛顿-拉夫逊潮流计算程序 (极坐标形式)
Newton-Raphson Power Flow Calculation Program

作者: 电力系统分析
日期: 2025-12-18
"""

import numpy as np
import time
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


class PowerFlowNR:
    """牛顿-拉夫逊潮流计算类"""
    
    def __init__(self, case_data):
        """
        初始化潮流计算
        
        参数:
            case_data: 包含bus, gen, branch, baseMVA的字典
        """
        self.baseMVA = case_data['baseMVA']
        self.bus = np.array(case_data['bus'], dtype=float)
        self.gen = np.array(case_data['gen'], dtype=float)
        self.branch = np.array(case_data['branch'], dtype=float)
        
        self.n_bus = len(self.bus)
        self.n_branch = len(self.branch)
        self.n_gen = len(self.gen)
        
        # 节点类型: 1-PQ, 2-PV, 3-Slack
        self.bus_type = self.bus[:, 1].astype(int)
        
        # 初始化状态变量
        self.V = self.bus[:, 7].copy()  # 电压幅值 (p.u.)
        self.theta = self.bus[:, 8].copy() * np.pi / 180  # 电压相角 (rad)
        
        # 导纳矩阵
        self.Ybus = None
        self.G = None  # 电导矩阵
        self.B = None  # 电纳矩阵
        
        # 节点注入功率
        self.Pgen = np.zeros(self.n_bus)
        self.Qgen = np.zeros(self.n_bus)
        self.Pload = self.bus[:, 2] / self.baseMVA
        self.Qload = self.bus[:, 3] / self.baseMVA
        
        # 收敛判据
        self.tol = 1e-8  # 收敛精度
        self.max_iter = 100  # 最大迭代次数
        
        # 计算结果
        self.converged = False
        self.iterations = 0
        self.runtime = 0
        
    def make_Ybus(self):
        """构建节点导纳矩阵"""
        Ybus = np.zeros((self.n_bus, self.n_bus), dtype=complex)
        
        for k in range(self.n_branch):
            fbus = int(self.branch[k, 0]) - 1  # From bus (转为0-based索引)
            tbus = int(self.branch[k, 1]) - 1  # To bus
            r = self.branch[k, 2]  # 电阻
            x = self.branch[k, 3]  # 电抗
            b = self.branch[k, 4]  # 充电电纳
            ratio = self.branch[k, 8]  # 变压器变比
            angle = self.branch[k, 9] * np.pi / 180  # 移相角
            status = self.branch[k, 10]  # 支路状态
            
            if status == 0:  # 支路断开
                continue
            
            # 串联导纳
            y = 1 / (r + 1j * x)
            
            # 并联充电导纳
            yc = 1j * b / 2
            
            if ratio == 0:
                ratio = 1.0
            
            # 变压器变比的复数形式
            tap = ratio * np.exp(1j * angle)
            
            # 构建导纳矩阵
            Ybus[fbus, fbus] += y / (tap * np.conj(tap)) + yc
            Ybus[tbus, tbus] += y + yc
            Ybus[fbus, tbus] += -y / np.conj(tap)
            Ybus[tbus, fbus] += -y / tap
        
        self.Ybus = Ybus
        self.G = Ybus.real
        self.B = Ybus.imag
        
    def setup_gen(self):
        """设置发电机注入功率"""
        for k in range(self.n_gen):
            bus_idx = int(self.gen[k, 0]) - 1
            if self.gen[k, 7] == 1:  # 发电机在线
                self.Pgen[bus_idx] += self.gen[k, 1] / self.baseMVA
                self.Qgen[bus_idx] += self.gen[k, 2] / self.baseMVA
                # PV节点和Slack节点设置电压幅值
                if self.bus_type[bus_idx] in [2, 3]:
                    self.V[bus_idx] = self.gen[k, 5]
    
    def calc_power(self):
        """计算节点注入功率"""
        P = np.zeros(self.n_bus)
        Q = np.zeros(self.n_bus)
        
        for i in range(self.n_bus):
            for j in range(self.n_bus):
                theta_ij = self.theta[i] - self.theta[j]
                P[i] += self.V[i] * self.V[j] * (
                    self.G[i, j] * np.cos(theta_ij) + 
                    self.B[i, j] * np.sin(theta_ij)
                )
                Q[i] += self.V[i] * self.V[j] * (
                    self.G[i, j] * np.sin(theta_ij) - 
                    self.B[i, j] * np.cos(theta_ij)
                )
        
        return P, Q
    
    def calc_jacobian(self):
        """计算雅可比矩阵"""
        # 分离PQ节点和PV节点
        pq_buses = np.where(self.bus_type == 1)[0]
        pv_buses = np.where(self.bus_type == 2)[0]
        
        n_pq = len(pq_buses)
        n_pv = len(pv_buses)
        n_total = n_pq + n_pv
        
        # 雅可比矩阵维度: (2*n_pq + n_pv) x (2*n_pq + n_pv)
        J = np.zeros((n_total + n_pq, n_total + n_pq))
        
        # 组合PQ和PV节点索引
        pqpv_buses = np.concatenate([pq_buses, pv_buses])
        
        # J11: ∂P/∂θ
        for i_idx, i in enumerate(pqpv_buses):
            for j_idx, j in enumerate(pqpv_buses):
                if i == j:
                    # 对角元素
                    val = 0
                    for k in range(self.n_bus):
                        if k != i:
                            theta_ik = self.theta[i] - self.theta[k]
                            val += self.V[i] * self.V[k] * (
                                self.G[i, k] * np.sin(theta_ik) - 
                                self.B[i, k] * np.cos(theta_ik)
                            )
                    J[i_idx, j_idx] = -val
                else:
                    # 非对角元素
                    theta_ij = self.theta[i] - self.theta[j]
                    J[i_idx, j_idx] = self.V[i] * self.V[j] * (
                        self.G[i, j] * np.sin(theta_ij) - 
                        self.B[i, j] * np.cos(theta_ij)
                    )
        
        # J12: ∂P/∂V (仅对PQ节点)
        for i_idx, i in enumerate(pqpv_buses):
            for j_idx, j in enumerate(pq_buses):
                if i == j:
                    # 对角元素
                    val = 0
                    for k in range(self.n_bus):
                        theta_ik = self.theta[i] - self.theta[k]
                        val += self.V[k] * (
                            self.G[i, k] * np.cos(theta_ik) + 
                            self.B[i, k] * np.sin(theta_ik)
                        )
                    J[i_idx, n_total + j_idx] = val + self.V[i] * self.G[i, i]
                else:
                    # 非对角元素
                    theta_ij = self.theta[i] - self.theta[j]
                    J[i_idx, n_total + j_idx] = self.V[i] * (
                        self.G[i, j] * np.cos(theta_ij) + 
                        self.B[i, j] * np.sin(theta_ij)
                    )
        
        # J21: ∂Q/∂θ (仅对PQ节点)
        for i_idx, i in enumerate(pq_buses):
            for j_idx, j in enumerate(pqpv_buses):
                if i == j:
                    # 对角元素
                    val = 0
                    for k in range(self.n_bus):
                        if k != i:
                            theta_ik = self.theta[i] - self.theta[k]
                            val += self.V[i] * self.V[k] * (
                                self.G[i, k] * np.cos(theta_ik) + 
                                self.B[i, k] * np.sin(theta_ik)
                            )
                    J[n_total + i_idx, j_idx] = val
                else:
                    # 非对角元素
                    theta_ij = self.theta[i] - self.theta[j]
                    J[n_total + i_idx, j_idx] = -self.V[i] * self.V[j] * (
                        self.G[i, j] * np.cos(theta_ij) + 
                        self.B[i, j] * np.sin(theta_ij)
                    )
        
        # J22: ∂Q/∂V (仅对PQ节点)
        for i_idx, i in enumerate(pq_buses):
            for j_idx, j in enumerate(pq_buses):
                if i == j:
                    # 对角元素
                    val = 0
                    for k in range(self.n_bus):
                        theta_ik = self.theta[i] - self.theta[k]
                        val += self.V[k] * (
                            self.G[i, k] * np.sin(theta_ik) - 
                            self.B[i, k] * np.cos(theta_ik)
                        )
                    J[n_total + i_idx, n_total + j_idx] = val - self.V[i] * self.B[i, i]
                else:
                    # 非对角元素
                    theta_ij = self.theta[i] - self.theta[j]
                    J[n_total + i_idx, n_total + j_idx] = self.V[i] * (
                        self.G[i, j] * np.sin(theta_ij) - 
                        self.B[i, j] * np.cos(theta_ij)
                    )
        
        return J, pq_buses, pv_buses
    
    def solve(self, verbose=True):
        """
        牛顿-拉夫逊法求解潮流
        
        参数:
            verbose: 是否打印迭代过程
        
        返回:
            converged: 是否收敛
            iterations: 迭代次数
        """
        start_time = time.time()
        
        # 构建导纳矩阵
        self.make_Ybus()
        
        # 设置发电机
        self.setup_gen()
        
        # 节点注入功率设定值
        Pspec = self.Pgen - self.Pload
        Qspec = self.Qgen - self.Qload
        
        # 分离PQ节点和PV节点
        pq_buses = np.where(self.bus_type == 1)[0]
        pv_buses = np.where(self.bus_type == 2)[0]
        pqpv_buses = np.concatenate([pq_buses, pv_buses])
        
        if verbose:
            print("\n" + "="*80)
            print("开始牛顿-拉夫逊潮流计算".center(70))
            print("="*80)
            print(f"系统规模: {self.n_bus} 节点, {self.n_branch} 支路, {self.n_gen} 发电机")
            print(f"PQ节点: {len(pq_buses)}, PV节点: {len(pv_buses)}, Slack节点: 1")
            print(f"收敛精度: {self.tol}, 最大迭代次数: {self.max_iter}")
            print("="*80)
            print(f"{'迭代次数':<10} {'最大不平衡量(P)':<20} {'最大不平衡量(Q)':<20}")
            print("-"*80)
        
        # 迭代求解
        for iteration in range(self.max_iter):
            # 计算节点注入功率
            P_calc, Q_calc = self.calc_power()
            
            # 功率不平衡量
            dP = Pspec[pqpv_buses] - P_calc[pqpv_buses]
            dQ = Qspec[pq_buses] - Q_calc[pq_buses]
            
            # 组合不平衡量向量
            dF = np.concatenate([dP, dQ])
            
            # 检查收敛性
            max_dP = np.max(np.abs(dP)) if len(dP) > 0 else 0
            max_dQ = np.max(np.abs(dQ)) if len(dQ) > 0 else 0
            
            if verbose:
                print(f"{iteration+1:<10} {max_dP:<20.6e} {max_dQ:<20.6e}")
            
            if max_dP < self.tol and max_dQ < self.tol:
                self.converged = True
                self.iterations = iteration + 1
                break
            
            # 计算雅可比矩阵
            J, _, _ = self.calc_jacobian()
            
            # 求解修正方程: J * dx = dF
            dx = np.linalg.solve(J, dF)
            
            # 更新状态变量
            n_pqpv = len(pqpv_buses)
            n_pq = len(pq_buses)
            
            # 更新相角
            self.theta[pqpv_buses] += dx[:n_pqpv]
            
            # 更新电压幅值(仅PQ节点)
            self.V[pq_buses] += dx[n_pqpv:]
        
        self.runtime = time.time() - start_time
        
        if verbose:
            print("="*80)
            if self.converged:
                print(f"✓ 潮流计算收敛! 迭代次数: {self.iterations}, 计算时间: {self.runtime:.4f}秒")
            else:
                print(f"✗ 潮流计算未收敛! 达到最大迭代次数: {self.max_iter}")
            print("="*80)
        
        return self.converged, self.iterations
    
    def calc_line_flow(self):
        """计算支路潮流"""
        line_flow = []
        
        for k in range(self.n_branch):
            fbus = int(self.branch[k, 0]) - 1
            tbus = int(self.branch[k, 1]) - 1
            r = self.branch[k, 2]
            x = self.branch[k, 3]
            b = self.branch[k, 4]
            ratio = self.branch[k, 8]
            angle = self.branch[k, 9] * np.pi / 180
            status = self.branch[k, 10]
            
            if status == 0:
                line_flow.append([0, 0, 0, 0, 0, 0])
                continue
            
            if ratio == 0:
                ratio = 1.0
            
            tap = ratio * np.exp(1j * angle)
            y = 1 / (r + 1j * x)
            yc = 1j * b / 2
            
            # From端电压和相角
            Vf = self.V[fbus] * np.exp(1j * self.theta[fbus])
            Vt = self.V[tbus] * np.exp(1j * self.theta[tbus])
            
            # From端注入功率
            If = (Vf / tap - Vt) * y + Vf / tap * yc
            Sf = Vf * np.conj(If)
            
            # To端注入功率
            It = (Vt - Vf / np.conj(tap)) * y + Vt * yc
            St = Vt * np.conj(It)
            
            # 线路损耗（仅串联阻抗上的I²Z损耗，不含充电功率）
            # 方法：计算流过串联阻抗的电流
            I_series = (Vf / tap - Vt) * y
            Sloss = I_series * np.conj(I_series) * (r + 1j * x)
            
            line_flow.append([
                Sf.real * self.baseMVA,
                Sf.imag * self.baseMVA,
                St.real * self.baseMVA,
                St.imag * self.baseMVA,
                Sloss.real * self.baseMVA,
                Sloss.imag * self.baseMVA
            ])
        
        return np.array(line_flow)
    
    def get_results(self):
        """获取计算结果"""
        # 重新计算节点注入功率
        P_calc, Q_calc = self.calc_power()
        
        # 更新发电机有功和无功功率
        for k in range(self.n_gen):
            bus_idx = int(self.gen[k, 0]) - 1
            if self.gen[k, 7] == 1:  # 在线发电机
                if self.bus_type[bus_idx] == 3:  # Slack节点
                    # Slack节点同时更新P和Q
                    self.gen[k, 1] = (P_calc[bus_idx] + self.Pload[bus_idx]) * self.baseMVA
                    self.gen[k, 2] = (Q_calc[bus_idx] + self.Qload[bus_idx]) * self.baseMVA
                elif self.bus_type[bus_idx] == 2:  # PV节点
                    # PV节点仅更新Q
                    self.gen[k, 2] = (Q_calc[bus_idx] + self.Qload[bus_idx]) * self.baseMVA
        
        # 计算支路潮流
        line_flow = self.calc_line_flow()
        
        results = {
            'converged': self.converged,
            'iterations': self.iterations,
            'runtime': self.runtime,
            'V': self.V.copy(),
            'theta': self.theta.copy() * 180 / np.pi,  # 转换为度
            'P': P_calc * self.baseMVA,
            'Q': Q_calc * self.baseMVA,
            'line_flow': line_flow,
            'bus': self.bus,
            'gen': self.gen,
            'branch': self.branch,
            'baseMVA': self.baseMVA
        }
        
        return results
