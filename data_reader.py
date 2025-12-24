"""
MATPOWER数据文件读取模块
支持读取.m格式的MATPOWER算例文件
"""

import re
import numpy as np


def parse_matpower_case(filename):
    """
    解析MATPOWER的.m格式算例文件
    
    参数:
        filename: 算例文件路径
    
    返回:
        case_data: 包含baseMVA, bus, gen, branch的字典
    """
    
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 初始化数据结构
    case_data = {
        'baseMVA': 100,  # 默认值
        'bus': [],
        'gen': [],
        'branch': []
    }
    
    # 提取baseMVA
    match = re.search(r'mpc\.baseMVA\s*=\s*(\d+)', content)
    if match:
        case_data['baseMVA'] = float(match.group(1))
    
    # 提取bus数据
    bus_match = re.search(r'mpc\.bus\s*=\s*\[(.*?)\];', content, re.DOTALL)
    if bus_match:
        bus_str = bus_match.group(1)
        case_data['bus'] = parse_matrix(bus_str)
    
    # 提取gen数据
    gen_match = re.search(r'mpc\.gen\s*=\s*\[(.*?)\];', content, re.DOTALL)
    if gen_match:
        gen_str = gen_match.group(1)
        case_data['gen'] = parse_matrix(gen_str)
    
    # 提取branch数据
    branch_match = re.search(r'mpc\.branch\s*=\s*\[(.*?)\];', content, re.DOTALL)
    if branch_match:
        branch_str = branch_match.group(1)
        case_data['branch'] = parse_matrix(branch_str)
    
    return case_data


def parse_matrix(matrix_str):
    """
    解析矩阵字符串为numpy数组
    
    参数:
        matrix_str: 矩阵字符串
    
    返回:
        numpy数组
    """
    # 移除注释
    lines = matrix_str.split('\n')
    data_lines = []
    
    for line in lines:
        # 移除注释
        line = re.sub(r'%.*', '', line)
        line = line.strip()
        
        if not line:
            continue
        
        # 移除前后的分号
        line = line.rstrip(';')
        
        # 分割数据
        values = line.split()
        if values:
            data_lines.append([float(v) for v in values])
    
    return data_lines


def load_case(case_name):
    """
    加载MATPOWER算例
    
    参数:
        case_name: 算例名称，如 'case5', 'case30', 'case118'
                  或完整文件路径
    
    返回:
        case_data: 算例数据字典
    """
    import os
    
    # 检查是否为完整路径
    if os.path.isfile(case_name):
        filename = case_name
    else:
        # 在潮流程序模板目录中查找
        template_dir = os.path.join(os.path.dirname(__file__), '潮流程序模板')
        filename = os.path.join(template_dir, f'{case_name}.m')
        
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"找不到算例文件: {filename}")
    
    print(f"正在加载算例: {filename}")
    case_data = parse_matpower_case(filename)
    
    n_bus = len(case_data['bus'])
    n_gen = len(case_data['gen'])
    n_branch = len(case_data['branch'])
    
    print(f"算例信息: {n_bus} 节点, {n_gen} 发电机, {n_branch} 支路")
    print(f"基准容量: {case_data['baseMVA']} MVA")
    
    return case_data


if __name__ == '__main__':
    # 测试数据读取
    import sys
    
    if len(sys.argv) > 1:
        case_name = sys.argv[1]
    else:
        case_name = 'case4gs'
    
    try:
        case_data = load_case(case_name)
        print("\n算例加载成功!")
        print(f"Bus数据样例 (前3行):")
        for row in case_data['bus'][:3]:
            print(row)
    except Exception as e:
        print(f"错误: {e}")
