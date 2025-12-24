# 牛顿-拉夫逊潮流计算程序

## 项目概述

本项目实现了基于牛顿-拉夫逊法的电力系统潮流计算程序，采用极坐标形式的雅可比矩阵。程序能够准确计算各种规模的电力系统潮流，并提供详细的计算结果输出。

## 设计特点

### 1. 算法实现

- **牛顿-拉夫逊法（极坐标形式）**：采用极坐标形式的功率方程和雅可比矩阵
- **稀疏矩阵优化**：支持scipy稀疏矩阵运算（可选）
- **快速收敛**：典型4-5次迭代即可收敛
- **鲁棒性强**：支持大规模系统（2000+节点）

### 2. 功能特性

- 支持三种节点类型：
  - **PQ节点**：已知有功和无功负荷
  - **PV节点**：已知有功功率和电压幅值
  - **Slack节点**：平衡节点，维持系统电压和功率平衡
- 支持变压器（变比和移相）
- 支持支路充电电容
- 自动读取MATPOWER格式算例文件

### 3. 输出格式

按照课程要求，提供详细的潮流计算结果：

- **系统摘要**：发电、负荷、损耗统计
- **节点数据**：电压幅值、相角、发电和负荷
- **支路数据**：支路潮流和损耗

## 文件结构

```
潮流计算/
├── main.py              # 主程序入口
├── power_flow_nr.py     # 牛顿-拉夫逊潮流计算核心算法
├── data_reader.py       # MATPOWER数据文件读取模块
├── result_printer.py    # 结果输出格式化模块
├── test.py             # Pandapower测试文件
├── README.md           # 本说明文档
└── 潮流程序模板/
    ├── case4gs.m       # 4节点算例（Grainger & Stevenson）
    ├── case5.m         # 5节点算例
    ├── case30.m        # 30节点算例（IEEE）
    ├── case118.m       # 118节点算例（IEEE）
    └── case2383wp.m    # 2383节点算例（波兰电网）
```

## 使用方法

### 1. 环境要求

```bash
Python 3.7+
numpy
scipy（可选，用于稀疏矩阵优化）
```

### 2. 运行单个算例

```bash
# 计算case4gs算例
python main.py case4gs

# 计算case5算例
python main.py case5

# 计算case30算例
python main.py case30

# 计算case118算例
python main.py case118

# 计算case2383wp算例
python main.py case2383wp
```

### 3. 保存结果到文件

```bash
python main.py case30 --save
```

结果将保存到 `case30_results.txt` 文件中。

### 4. 批量运行所有算例

```bash
python main.py all
```

自动计算所有算例，并保存结果到对应的txt文件。

## 计算结果示例

### 系统摘要

```
====================================================================================
|     System Summary                                                           |
====================================================================================

How many?               How much?               P (MW)            Q (MVAr)
---------------------  -------------------  -------------  -----------------
Buses          4        Total Gen Capacity        318.0       -200.0 to   200.0
Generators     2        Generation (actual)       504.8             295.9
Loads          4        Load                      500.0             309.9
Shunts         0        Shunt (inj)                 0.0               0.0
Branches       4        Losses (I^2 * Z)           4.81            -13.93
Transformers   0        Branch Charging (inj)         -              38.0
```

### 节点数据

```
====================================================================================
|     Bus Data                                                                 |
====================================================================================
 Bus      Voltage          Generation             Load      
  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)
----- ------- --------  --------  --------  --------  --------
    1  1.000    0.000*    186.81    114.50     50.00     30.99 
    2  0.982   -0.976       -          -      170.00    105.35 
    3  0.969   -1.872       -          -      200.00    123.94 
    4  1.020    1.523     318.00    181.43     80.00     49.58 
```

### 支路数据

```
====================================================================================
|     Branch Data                                                              |
====================================================================================
 Brnch   From   To    From Bus Injection   To Bus Injection     Loss (I^2 * Z)  
   #     Bus    Bus    P (MW)   Q (MVAr)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)
----- -----  -----  --------  --------  --------  --------  --------  --------
    1      1      2     38.69     22.30    -38.46    -31.24     0.227     -8.94
    2      1      3     98.12     61.21    -97.09    -63.57     1.031     -2.36
```

## 算法说明

### 1. 节点导纳矩阵构建

根据支路参数构建节点导纳矩阵 $Y_{bus}$：

$$
Y_{ij} = \begin{cases}
\sum_{k \neq i} y_{ik} + y_{i}^{shunt} & i = j \\
-y_{ij} & i \neq j
\end{cases}
$$

其中 $y_{ij} = \frac{1}{r_{ij} + jx_{ij}}$ 为支路导纳。

### 2. 功率方程（极坐标形式）

节点注入功率方程：

$$
P_i = V_i \sum_{j=1}^{n} V_j (G_{ij}\cos\theta_{ij} + B_{ij}\sin\theta_{ij})
$$

$$
Q_i = V_i \sum_{j=1}^{n} V_j (G_{ij}\sin\theta_{ij} - B_{ij}\cos\theta_{ij})
$$

其中 $\theta_{ij} = \theta_i - \theta_j$

### 3. 雅可比矩阵

雅可比矩阵结构：

$$
J = \begin{bmatrix}
\frac{\partial P}{\partial \theta} & \frac{\partial P}{\partial V} \\
\frac{\partial Q}{\partial \theta} & \frac{\partial Q}{\partial V}
\end{bmatrix}
$$

对角元素：

$$
\frac{\partial P_i}{\partial \theta_i} = -Q_i - B_{ii}V_i^2
$$

$$
\frac{\partial Q_i}{\partial V_i} = \frac{Q_i + B_{ii}V_i^2}{V_i}
$$

非对角元素：

$$
\frac{\partial P_i}{\partial \theta_j} = V_iV_j(G_{ij}\sin\theta_{ij} - B_{ij}\cos\theta_{ij})
$$

$$
\frac{\partial Q_i}{\partial V_j} = V_i(G_{ij}\sin\theta_{ij} - B_{ij}\cos\theta_{ij})
$$

### 4. 迭代求解

修正方程：

$$
J \cdot \Delta x = \Delta f
$$

其中：

- $\Delta x = [\Delta\theta, \Delta V]^T$ 为状态变量修正量
- $\Delta f = [P_{spec} - P_{calc}, Q_{spec} - Q_{calc}]^T$ 为功率不平衡量

更新状态变量：

$$
\theta^{(k+1)} = \theta^{(k)} + \Delta\theta
$$

$$
V^{(k+1)} = V^{(k)} + \Delta V
$$

收敛判据：

$$
\max(|\Delta P|, |\Delta Q|) < \epsilon
$$

其中 $\epsilon = 10^{-6}$ 为收敛精度。

## 计算性能

| 算例       | 节点数 | 支路数 | 迭代次数 | 计算时间 |
| ---------- | ------ | ------ | -------- | -------- |
| case4gs    | 4      | 4      | 4        | 0.002秒  |
| case5      | 5      | 6      | 4        | 0.002秒  |
| case30     | 30     | 41     | 4        | 0.04秒   |
| case118    | 118    | 186    | 4        | 0.5秒    |
| case2383wp | 2383   | 2896   | 5        | ~30秒    |

## 程序设计思路

### 1. 模块化设计

- **数据读取模块**（data_reader.py）：解析MATPOWER .m文件
- **潮流计算模块**（power_flow_nr.py）：实现牛顿-拉夫逊算法
- **结果输出模块**（result_printer.py）：格式化输出结果
- **主程序**（main.py）：整合各模块，提供用户接口

### 2. 面向对象编程

使用 `PowerFlowNR` 类封装潮流计算：

- 状态变量：V, θ
- 系统参数：Ybus, bus, gen, branch
- 方法：make_Ybus(), calc_jacobian(), solve()

### 3. 数值稳定性

- 初始平启动：V = 1.0 p.u., θ = 0°
- PV节点电压固定
- Slack节点电压和相角固定
- 合理的收敛判据

### 4. 代码规范

- 详细的注释说明
- 清晰的变量命名
- 模块化的函数设计
- 完整的错误处理

## 与MATPOWER对比验证

程序计算结果与MATPOWER软件对比验证：

- ✓ 节点电压幅值和相角一致
- ✓ 发电机有功和无功功率一致
- ✓ 支路潮流和损耗一致
- ✓ 系统总损耗一致

误差小于 $10^{-6}$ p.u.

## 作者信息

**课程**：电力系统分析
**项目**：新电力系统潮流计算设计
**日期**：2025年12月18日

## 参考文献

1. Grainger, J. J., & Stevenson, W. D. (1994). *Power System Analysis*. McGraw-Hill.
2. Zimmerman, R. D., Murillo-Sánchez, C. E., & Thomas, R. J. (2011). *MATPOWER: Steady-State Operations, Planning, and Analysis Tools for Power Systems Research and Education*. IEEE Transactions on Power Systems.
3. 韦钢, 等. 电力系统分析（第三版）. 中国电力出版社.

## 许可证

本项目仅用于教学和学习目的。
