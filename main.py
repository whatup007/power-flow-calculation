"""
牛顿-拉夫逊潮流计算主程序
支持多种算例的潮流计算

运行方式:
    python main.py case4gs
    python main.py case5
    python main.py case30
    python main.py case118
    python main.py case2383wp
"""

import sys
import os
from data_reader import load_case
from power_flow_nr import PowerFlowNR
from result_printer import print_results, save_results_to_file


def main(case_name='case4gs', save_to_file=False, use_sparse=False, outdir=None):
    """
    主函数
    
    参数:
        case_name: 算例名称
        save_to_file: 是否保存结果到文件
    """
    print("\n" + "="*84)
    print("牛顿-拉夫逊潮流计算程序 (极坐标形式)".center(74))
    print("Newton-Raphson Power Flow Calculation Program".center(74))
    print("="*84)
    
    try:
        # 加载算例数据
        case_data = load_case(case_name)
        
        # 创建潮流计算对象
        pf = PowerFlowNR(case_data, use_sparse=use_sparse)
        
        # 求解潮流
        converged, iterations = pf.solve(verbose=True)
        
        # 获取结果
        results = pf.get_results()
        
        # 打印结果
        print_results(results)
        
        # 保存结果到文件
        # 需求：使用稀疏矩阵加速时，结果保存在当前目录，命名为 <case_name>_sparse_result.txt
        if use_sparse:
            output_file = f"{case_name}_sparse_result.txt"
            save_results_to_file(results, output_file)
        elif save_to_file:
            # 非稀疏：沿用原有保存逻辑，可选 outdir
            target_dir = outdir
            if target_dir:
                os.makedirs(target_dir, exist_ok=True)
                output_file = os.path.join(target_dir, f"{case_name}_results.txt")
            else:
                output_file = f"{case_name}_results.txt"
            save_results_to_file(results, output_file)
        
        return results
        
    except FileNotFoundError as e:
        print(f"\n错误: {e}")
        print("\n可用的算例:")
        print("  - case4gs")
        print("  - case5")
        print("  - case30")
        print("  - case118")
        print("  - case2383wp")
        return None
    except Exception as e:
        print(f"\n计算过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        return None


def run_all_cases(use_sparse=False):
    """运行所有算例
    
    参数:
        use_sparse: 是否使用稀疏矩阵
    """
    cases = ['case4gs', 'case5', 'case30', 'case118', 'case2383wp']
    
    print("\n" + "="*84)
    if use_sparse:
        print("批量运行所有算例 (稀疏矩阵加速)".center(74))
    else:
        print("批量运行所有算例".center(74))
    print("="*84)
    
    results_summary = []
    
    for case_name in cases:
        print(f"\n\n{'='*84}")
        print(f"正在计算算例: {case_name}".center(74))
        print("="*84)
        
        try:
            # 加载算例
            case_data = load_case(case_name)
            
            # 求解潮流
            pf = PowerFlowNR(case_data, use_sparse=use_sparse)
            converged, iterations = pf.solve(verbose=False)
            
            results_summary.append({
                'case': case_name,
                'converged': converged,
                'iterations': iterations,
                'runtime': pf.runtime,
                'n_bus': pf.n_bus,
                'n_branch': pf.n_branch
            })
            
            if converged:
                print(f"✓ {case_name} 计算成功: 迭代{iterations}次, 耗时{pf.runtime:.4f}秒")
            else:
                print(f"✗ {case_name} 未收敛")
            
            # 保存结果
            results = pf.get_results()
            if use_sparse:
                output_file = f"{case_name}_sparse_result.txt"
            else:
                output_file = f"{case_name}_results.txt"
            save_results_to_file(results, output_file)
            
        except Exception as e:
            print(f"✗ {case_name} 计算失败: {e}")
            results_summary.append({
                'case': case_name,
                'converged': False,
                'iterations': 0,
                'runtime': 0,
                'n_bus': 0,
                'n_branch': 0
            })
    
    # 打印汇总信息
    print("\n\n" + "="*84)
    print("所有算例计算汇总".center(74))
    print("="*84)
    print(f"{'算例':<15} {'节点数':<10} {'支路数':<10} {'收敛':<10} {'迭代次数':<12} {'计算时间(秒)':<15}")
    print("-"*84)
    
    for r in results_summary:
        status = "✓" if r['converged'] else "✗"
        print(f"{r['case']:<15} {r['n_bus']:<10} {r['n_branch']:<10} "
              f"{status:<10} {r['iterations']:<12} {r['runtime']:<15.4f}")
    
    print("="*84)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == 'all':
            # 运行所有算例
            args = sys.argv[2:]
            use_sparse = ('--sparse' in args)
            run_all_cases(use_sparse=use_sparse)
        else:
            # 运行指定算例并解析可选参数
            case_name = sys.argv[1]
            args = sys.argv[2:]
            save_file = ('--save' in args) or ('-s' in args)
            use_sparse = ('--sparse' in args)
            # 解析输出目录
            outdir = None
            if '--outdir' in args:
                idx = args.index('--outdir')
                if idx + 1 < len(args):
                    outdir = args[idx + 1]
            main(case_name, save_to_file=save_file, use_sparse=use_sparse, outdir=outdir)
    else:
        # 默认运行case4gs算例
        print("\n使用方法:")
        print("  python main.py <case_name>                           # 计算指定算例")
        print("  python main.py <case_name> --sparse                  # 稀疏加速，并自动保存到当前目录为 <case_name>_sparse_result.txt")
        print("  python main.py <case_name> --save                    # 非稀疏下保存结果到当前目录")
        print("  python main.py <case_name> --save --outdir <path>    # 非稀疏下保存到指定目录")
        print("  python main.py all                                    # 计算所有算例")
        print("  python main.py all --sparse                          # 计算所有算例（使用稀疏矩阵加速）")
        print("\n可用算例: case4gs, case5, case30, case118, case2383wp")
        print("\n现在运行默认算例: case4gs\n")
        main('case4gs')
