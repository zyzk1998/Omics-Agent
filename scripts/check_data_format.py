#!/usr/bin/env python3
"""
检查数据文件格式
用于诊断工具执行失败的原因
"""
import sys
import os
import pandas as pd

def check_data_format(file_path):
    """检查数据文件格式"""
    print("=" * 80)
    print("🔍 步骤4：检查数据文件格式")
    print("=" * 80)
    
    if not os.path.exists(file_path):
        print(f"❌ 文件不存在: {file_path}")
        return False
    
    try:
        print(f"\n1. 读取数据文件: {file_path}")
        df = pd.read_csv(file_path, nrows=10)
        print(f"   ✅ 数据文件读取成功")
        print(f"   - 形状: {df.shape}")
        print(f"   - 列名: {df.columns.tolist()}")
        
        print(f"\n2. 列类型分析:")
        numeric_cols = []
        non_numeric_cols = []
        for col in df.columns:
            dtype = df[col].dtype
            n_unique = df[col].nunique()
            col_type = "数值型" if pd.api.types.is_numeric_dtype(df[col]) else "非数值型"
            
            if pd.api.types.is_numeric_dtype(df[col]):
                numeric_cols.append((col, dtype, n_unique))
            else:
                non_numeric_cols.append((col, dtype, n_unique))
            
            print(f"   - {col}: {dtype} ({col_type}, 唯一值: {n_unique})")
        
        print(f"\n3. 分组列检测:")
        if non_numeric_cols:
            print(f"   ✅ 找到 {len(non_numeric_cols)} 个非数值列:")
            for col, dtype, n_unique in non_numeric_cols:
                print(f"      - {col}: {dtype} (唯一值: {n_unique})")
        else:
            print(f"   ⚠️ 未找到非数值列（所有列都是数值型）")
        
        # 检查可能的分组列（唯一值在2-10之间）
        potential_group_cols = []
        for col in df.columns:
            n_unique = df[col].nunique()
            if 2 <= n_unique <= 10:
                potential_group_cols.append((col, n_unique, df[col].dtype))
        
        if potential_group_cols:
            print(f"\n   ✅ 找到 {len(potential_group_cols)} 个可能的分组列（唯一值在2-10之间）:")
            for col, n_unique, dtype in potential_group_cols:
                unique_values = df[col].unique()[:5].tolist()
                print(f"      - {col}: {dtype} (唯一值: {n_unique}, 示例: {unique_values})")
        else:
            print(f"\n   ⚠️ 未找到可能的分组列（唯一值在2-10之间）")
        
        # 检查特定的分组列名
        common_group_names = ['Sample', 'Group', 'Condition', 'Treatment', 'Class', 'Label']
        found_group_cols = []
        for col in df.columns:
            if col in common_group_names:
                found_group_cols.append(col)
        
        if found_group_cols:
            print(f"\n   ✅ 找到常见分组列名:")
            for col in found_group_cols:
                n_unique = df[col].nunique()
                unique_values = df[col].unique()[:5].tolist()
                print(f"      - {col}: 唯一值: {n_unique}, 示例: {unique_values}")
        else:
            print(f"\n   ⚠️ 未找到常见分组列名（Sample, Group, Condition, Treatment, Class, Label）")
        
        print(f"\n4. 数据格式评估:")
        if non_numeric_cols:
            print(f"   ✅ 数据格式正常（包含非数值列）")
            if potential_group_cols:
                print(f"   ✅ 找到可能的分组列")
            else:
                print(f"   ⚠️ 非数值列的唯一值不在2-10之间，可能不适合作为分组列")
        else:
            print(f"   ❌ 数据格式问题：所有列都是数值型，缺少分组信息列")
            print(f"   💡 建议：")
            print(f"      1. 在数据文件中添加一列分组信息（如 'Group', 'Condition' 等）")
            print(f"      2. 或者使用其他列作为分组列（如果唯一值在2-10之间）")
        
        return True
        
    except Exception as e:
        print(f"❌ 数据文件读取失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """主函数"""
    # 默认文件路径
    default_path = "guest/20260128_112036/cow_diet.csv"
    
    if len(sys.argv) > 1:
        file_path = sys.argv[1]
    else:
        file_path = default_path
        print(f"⚠️ 未指定文件路径，使用默认路径: {file_path}")
        print(f"   使用方法: python {sys.argv[0]} <文件路径>")
        print()
    
    check_data_format(file_path)

if __name__ == "__main__":
    main()
