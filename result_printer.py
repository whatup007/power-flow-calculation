"""
潮流计算结果输出模块
按照课程要求格式输出潮流计算结果
"""

import numpy as np


def print_system_summary(results):
    """打印系统摘要信息"""
    bus = results['bus']
    gen = results['gen']
    branch = results['branch']
    baseMVA = results['baseMVA']
    line_flow = results['line_flow']
    V = results['V']
    theta = results['theta']
    
    n_bus = len(bus)
    n_gen = len(gen)
    n_branch = len(branch)
    
    # 统计发电机容量
    gen_p_total = 0
    gen_q_min = 0
    gen_q_max = 0
    gen_p_actual = 0
    gen_q_actual = 0
    
    for g in gen:
        if g[7] == 1:  # 在线发电机
            gen_p_total += g[8]  # Pmax
            gen_q_min += g[4]    # Qmin
            gen_q_max += g[3]    # Qmax
            gen_p_actual += g[1] # Pg
            gen_q_actual += g[2] # Qg
    
    # 统计负荷
    load_p = np.sum(bus[:, 2])
    load_q = np.sum(bus[:, 3])
    
    # 统计并联电容器/电抗器
    shunt_p = np.sum(bus[:, 4]) * baseMVA / 1.0  # Gs * V^2
    shunt_q = np.sum(bus[:, 5]) * baseMVA / 1.0  # Bs * V^2
    
    # 计算损耗
    loss_p = np.sum(line_flow[:, 4])  # I²R有功损耗
    loss_q = np.sum(line_flow[:, 5])  # I²X无功损耗
    
    # 支路充电无功
    branch_charging = 0
    for k in range(n_branch):
        if branch[k, 10] == 1:  # 在线支路
            b = branch[k, 4]
            fbus = int(branch[k, 0]) - 1
            tbus = int(branch[k, 1]) - 1
            Vf = V[fbus]
            Vt = V[tbus]
            branch_charging += b * (Vf**2 + Vt**2) / 2 * baseMVA
    
    # 统计变压器
    n_transformer = 0
    for k in range(n_branch):
        if branch[k, 8] != 0:  # ratio != 0表示变压器
            n_transformer += 1
    
    # 电压幅值统计
    v_min = np.min(V)
    v_max = np.max(V)
    v_min_bus = np.argmin(V) + 1
    v_max_bus = np.argmax(V) + 1
    
    # 电压相角统计
    theta_min = np.min(theta)
    theta_max = np.max(theta)
    theta_min_bus = np.argmin(theta) + 1
    theta_max_bus = np.argmax(theta) + 1
    
    # 支路损耗统计
    p_loss_max = np.max(line_flow[:, 4]) if n_branch > 0 else 0
    q_loss_max = np.max(line_flow[:, 5]) if n_branch > 0 else 0
    p_loss_max_line = np.argmax(line_flow[:, 4]) + 1 if n_branch > 0 else 0
    q_loss_max_line = np.argmax(line_flow[:, 5]) + 1 if n_branch > 0 else 0
    
    print("\n" + "="*84)
    print("|     System Summary                                                           |")
    print("="*84)
    print()
    print(f"{'How many?':<24}{'How much?':<24}{'P (MW)':<18}Q (MVAr)")
    print(f"{'-'*21:21}  {'-'*19:19}  {'-'*13:13}  {'-'*17:17}")
    print(f"{'Buses':<15}{n_bus:<9}{'Total Gen Capacity':<23}{gen_p_total:>8.1f}      {gen_q_min:>7.1f} to {gen_q_max:>7.1f}")
    print(f"{'Generators':<15}{n_gen:<9}{'Generation (actual)':<23}{gen_p_actual:>8.1f}      {gen_q_actual:>12.1f}")
    print(f"{'Loads':<15}{np.sum(bus[:, 2] > 0):<9.0f}{'Load':<23}{load_p:>8.1f}      {load_q:>12.1f}")
    print(f"{'Shunts':<15}{np.sum((bus[:, 4] != 0) | (bus[:, 5] != 0)):<9.0f}{'Shunt (inj)':<23}{shunt_p:>8.1f}      {shunt_q:>12.1f}")
    print(f"{'Branches':<15}{n_branch:<9}{'Losses (I^2 * R)':<23}{loss_p:>8.2f}      {'':>12}")
    print(f"{'':15}{'':9}{'Losses (I^2 * X)':<23}{'     -':>8}      {loss_q:>12.2f}")
    print(f"{'Transformers':<15}{n_transformer:<9}{'Branch Charging (inj)':<23}{'     -':>8}      {branch_charging:>12.1f}")
    print()
    print()
    print(f"{'':26}Minimum                      Maximum")
    print(f"{' '*18}{'-'*25}  {'-'*32}")
    print(f"Voltage Magnitude   {v_min:5.3f} p.u. @ bus {v_min_bus:<9}{v_max:5.3f} p.u. @ bus {v_max_bus:<5}")
    print(f"Voltage Angle      {theta_min:6.2f} deg   @ bus {theta_min_bus:<9}{theta_max:5.2f} deg   @ bus {theta_max_bus:<5}")
    print(f"P Losses (I^2*R)             -                  {p_loss_max:>6.2f} MW    @ line {int(branch[p_loss_max_line-1, 0])}-{int(branch[p_loss_max_line-1, 1])}")
    print(f"Q Losses (I^2*X)             -                  {q_loss_max:>6.2f} MVAr  @ line {int(branch[q_loss_max_line-1, 0])}-{int(branch[q_loss_max_line-1, 1])}")
    print()


def print_bus_data(results):
    """打印节点数据"""
    bus = results['bus']
    gen = results['gen']
    V = results['V']
    theta = results['theta']
    P = results['P']
    Q = results['Q']
    baseMVA = results['baseMVA']
    
    n_bus = len(bus)
    
    # 计算每个节点的发电和负荷
    bus_gen_p = np.zeros(n_bus)
    bus_gen_q = np.zeros(n_bus)
    
    for g in gen:
        if g[7] == 1:  # 在线发电机
            bus_idx = int(g[0]) - 1
            bus_gen_p[bus_idx] += g[1]
            bus_gen_q[bus_idx] += g[2]
    
    bus_load_p = bus[:, 2]
    bus_load_q = bus[:, 3]
    
    print("="*84)
    print("|     Bus Data                                                                 |")
    print("="*84)
    print(" Bus      Voltage          Generation             Load        ")
    print("  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)")
    print("----- ------- --------  --------  --------  --------  --------")
    
    total_gen_p = 0
    total_gen_q = 0
    total_load_p = 0
    total_load_q = 0
    
    for i in range(n_bus):
        bus_num = int(bus[i, 0])
        bus_type = int(bus[i, 1])
        
        # 确定节点类型标记
        if bus_type == 3:  # Slack
            type_mark = '*'
        else:
            type_mark = ' '
        
        # 发电机输出 (检查是否有发电机或为PV/Slack节点)
        has_gen = (bus_gen_p[i] != 0) or (bus_type in [2, 3] and any(gen[:, 0] == bus_num))
        if has_gen:
            gen_p_str = f"{bus_gen_p[i]:>8.2f}"
            gen_q_str = f"{bus_gen_q[i]:>8.2f}"
        else:
            gen_p_str = "    -    "
            gen_q_str = "    -    "
        
        # 负荷
        if bus_load_p[i] > 0:
            load_p_str = f"{bus_load_p[i]:>8.2f}"
            load_q_str = f"{bus_load_q[i]:>8.2f}"
        else:
            load_p_str = "    -    "
            load_q_str = "    -    "
        
        print(f"{bus_num:>5}  {V[i]:>5.3f}  {theta[i]:>7.3f}{type_mark}  "
              f"{gen_p_str}  {gen_q_str}  {load_p_str}  {load_q_str}")
        
        total_gen_p += bus_gen_p[i]
        total_gen_q += bus_gen_q[i]
        total_load_p += bus_load_p[i]
        total_load_q += bus_load_q[i]
    
    print(" "*33 + "--------  --------  --------  --------")
    print(f"{'Total:':>33}  {total_gen_p:>8.2f}  {total_gen_q:>8.2f}  "
          f"{total_load_p:>8.2f}  {total_load_q:>8.2f}")
    print()


def print_branch_data(results):
    """打印支路数据"""
    branch = results['branch']
    line_flow = results['line_flow']
    
    n_branch = len(branch)
    
    print("="*84)
    print("|     Branch Data                                                              |")
    print("="*84)
    print(" Brnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  ")
    print("   #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)")
    print("----- -----  -----  --------  --------  --------  --------  --------  --------")
    
    total_loss_p = 0
    total_loss_q = 0
    
    for k in range(n_branch):
        if branch[k, 10] == 0:  # 支路断开
            continue
        
        fbus = int(branch[k, 0])
        tbus = int(branch[k, 1])
        
        pf = line_flow[k, 0]
        qf = line_flow[k, 1]
        pt = line_flow[k, 2]
        qt = line_flow[k, 3]
        ploss = line_flow[k, 4]
        qloss = line_flow[k, 5]
        
        print(f"{k+1:>5}  {fbus:>5}  {tbus:>5}  "
              f"{pf:>8.2f}  {qf:>8.2f}  {pt:>8.2f}  {qt:>8.2f}  "
              f"{ploss:>8.3f}  {qloss:>8.2f}")
        
        total_loss_p += ploss
        total_loss_q += qloss
    
    print(" "*62 + "--------  --------")
    print(f"{'Total:':>62}  {total_loss_p:>8.3f}  {total_loss_q:>8.2f}")
    print()


def print_results(results):
    """
    打印完整的潮流计算结果
    
    参数:
        results: 包含所有计算结果的字典
    """
    converged = results['converged']
    iterations = results['iterations']
    runtime = results['runtime']
    
    print("\n" + "="*84)
    print("|     NEWTON-RAPHSON POWER FLOW CALCULATION RESULTS                          |")
    print("="*84)
    
    if converged:
        print(f"\n✓ 潮流计算收敛")
        print(f"  迭代次数: {iterations}")
        print(f"  计算时间: {runtime:.4f} 秒")
    else:
        print(f"\n✗ 潮流计算未收敛")
        print(f"  已达到最大迭代次数")
        return
    
    # 打印系统摘要
    print_system_summary(results)
    
    # 打印节点数据
    print_bus_data(results)
    
    # 打印支路数据
    print_branch_data(results)
    
    print("="*84)


def save_results_to_file(results, filename):
    """
    将结果保存到文件
    
    参数:
        results: 计算结果字典
        filename: 输出文件名
    """
    import sys
    from io import StringIO
    
    # 重定向stdout到字符串
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    # 打印结果
    print_results(results)
    
    # 获取输出内容
    output = sys.stdout.getvalue()
    
    # 恢复stdout
    sys.stdout = old_stdout
    
    # 写入文件
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(output)
    
    print(f"\n结果已保存到文件: {filename}")


if __name__ == '__main__':
    # 测试输出格式
    print("输出格式测试")
    print("="*84)
    print("此模块用于格式化输出潮流计算结果")
    print("请运行 main.py 查看完整输出")
