import os
import glob
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from collections import Counter

# ==========================================
# 1. 核心参数与全局变量配置
# ==========================================
# 存放 10万个 SNP 文件的总目录
INPUT_DIR = "../ano/"
FILE_EXTENSION = "*.ano" 
SEPARATOR = "\t"

# 自定义的列名 (无表头文件必备)
COL_NAMES = [
    'genome_position', 'ref', 'alt', 'codon_position', 
    'ano', 'codon_change', 'locus_tag', 'gene', 'function', 'group'
]

# 在读取时只加载这 4 列，极大节省内存
USE_COLS = ['ref', 'alt', 'codon_change', 'gene']

# 耐药基因列表 (转换为 set 提升查询速度)
DR_GENES = {"rpoB", "katG", "inhA", "pncA", "gyrA", "gyrB", \
    "embB", "ethA", "Rv0678", "aptE", "pepQ", "rrs", "gid", "rpsL", "eis","tlyA"}

# FFD 密码子前两个碱基的组合
FFD_PREFIXES = {"CT", "GT", "TC", "CC", "AC", "GC", "CG", "GG"}

# ==========================================
# 2. 定义单个文件的处理函数 (Worker Function)
# ==========================================
def process_single_file(file_path):
    local_counts = {
        'A': Counter(), 'C': Counter(), 
        'G': Counter(), 'T': Counter()
    }
    
    try:
        # 读取文件：指定无表头 (header=None)，并赋予列名
        df = pd.read_csv(
            file_path, 
            sep=SEPARATOR, 
            header=None, 
            names=COL_NAMES, 
            usecols=USE_COLS, 
            dtype=str
        )
        
        # 1. 排除耐药基因 (处理 NaN 的情况)
        df = df[~df['gene'].isin(DR_GENES)]
        
        # 2. 仅保留单碱基替换 (SNPs)
        df = df[(df['ref'].str.len() == 1) & (df['alt'].str.len() == 1)]
        
        # 3. 过滤 codon_change 列为空的数据
        df = df.dropna(subset=['codon_change'])
        
        # 4. 【核心优化】严谨且高速的 FFD 过滤逻辑
        # 格式为 "CAG-CGG" (长度必须为 7)
        df = df[df['codon_change'].str.len() == 7]
        
        # 提取参考和突变密码子的前两位
        # C A G - C G G 
        # 0 1 2 3 4 5 6  (索引切片：[0:2]是参考前两位，[4:6]是突变前两位)
        ref_prefix = df['codon_change'].str[:2]
        alt_prefix = df['codon_change'].str[4:6]
        
        # FFD 条件 A：前两位必须属于那 8 种特殊的氨基酸密码子
        cond_ffd = ref_prefix.isin(FFD_PREFIXES)
        
        # FFD 条件 B：突变前后，密码子的前两位必须完全一致 (意味着突变只发生在第3位，绝对是同义突变)
        cond_synonymous = (ref_prefix == alt_prefix)
        
        # 应用过滤条件
        df = df[cond_ffd & cond_synonymous]
        
        # 5. 统计最终留下的真实 FFD 突变
        for r, a in zip(df['ref'], df['alt']):
            if r in local_counts and a in local_counts:
                local_counts[r][a] += 1
                
    except Exception as e:
        # 忽略损坏的文件，保证主进程继续运行
        pass
        
    return local_counts

# ==========================================
# 3. 主程序：多进程并发执行与结果汇总
# ==========================================
if __name__ == '__main__':
    print(f"1. 正在扫描目录 {INPUT_DIR} 寻找突变文件...")
    all_files = glob.glob(os.path.join(INPUT_DIR, f"**/{FILE_EXTENSION}"), recursive=True)
    total_files = len(all_files)
    print(f"   找到 {total_files} 个文件准备处理。\n")

    if total_files == 0:
        print("未找到任何文件，请检查 INPUT_DIR 和 FILE_EXTENSION。")
        exit()

    global_mutation_counts = {
        'A': {'C': 0, 'G': 0, 'T': 0},
        'C': {'A': 0, 'G': 0, 'T': 0},
        'G': {'A': 0, 'C': 0, 'T': 0},
        'T': {'A': 0, 'C': 0, 'G': 0}
    }

    # 预留 2 个核心，防止服务器卡死
    num_workers = max(1, cpu_count() - 2)
    print(f"2. 启动多进程并发处理 (开启 {num_workers} 个 CPU 核心)...")

    with Pool(processes=num_workers) as pool:
        for local_counts in tqdm(pool.imap_unordered(process_single_file, all_files), total=total_files, desc="Processing Samples"):
            for ref_base, alt_counts in local_counts.items():
                for alt_base, count in alt_counts.items():
                    global_mutation_counts[ref_base][alt_base] += count

    print("\n3. 所有样本处理完成！正在计算经验 GTR 概率矩阵...\n")

    # ==========================================
    # 4. 计算并输出最终结果
    # ==========================================
    print("# ================= GTR 经验概率矩阵 ================= #")
    for ref_base in ['A', 'C', 'G', 'T']:
        total_mutations_from_ref = sum(global_mutation_counts[ref_base].values())
        
        if total_mutations_from_ref == 0:
            print(f"# 警告: {ref_base} 碱基没有任何 FFD 突变记录！")
            continue
            
        for alt_base in ['A', 'C', 'G', 'T']:
            if ref_base != alt_base:
                count = global_mutation_counts[ref_base][alt_base]
                freq = count / total_mutations_from_ref
                print(f"# {ref_base} -> {alt_base}: count = {count:<10} frequency = {freq:.12f}")