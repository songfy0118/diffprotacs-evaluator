# DiffPROTACs Evaluator —— 给 DiffPROTACs 的评估系统

给 DiffPROTACs (扩散模型生成 PROTAC linker) 配的一套**评价指标 + 打分管线**。
不动 DiffPROTACs 的代码,作为独立工具放在旁边就能用。

> 解决的问题:DiffPROTACs 生成完只有 3D `.xyz` 文件,**没有 SMILES、没有 validity、没有 QED、没有 docking score**。这个工具把它们一次性全补齐。

---

## 这个工具做什么

```
DiffPROTACs 生成的 .xyz 文件
        │
        ▼  evaluator/xyz_to_smiles.py     (3D 坐标 → 化学键 → SMILES)
        │
        ▼  evaluator/scorers.py           (RDKit 打分:validity/QED/SA/Lipinski/alerts/novelty)
        │
        ▼  evaluator/docking_bridge.py    (可选:调 DockStream/Vina 做 docking)
        │
        ▼  CSV 报告 + Top-K + 终端 summary
```

输出每个分子的:
- **SMILES** (从 .xyz 重建出来的)
- **validity** —— 能否解析成合法分子
- **QED** —— 类药性 (Bickerton 2012)
- **SA proxy** —— 合成可及性代理
- **Lipinski fraction** —— Ro5 通过率 (PROTAC 一般违反,作为软信号)
- **size_window** —— 重原子数是否在合理范围
- **alerts_pass** —— 是否有 PAINS-like 不良基团 (硬门控)
- **novelty** —— 相对参考库(例如 PROTAC-DB 3.0)的新颖度
- **dock_score** —— Vina 打分 (可选,需要 DockStream)
- **total_score** —— 加权几何平均(与 REINVENT 的 `custom_product` 一致)

---

## 👋 给我朋友: 30 秒上手

```bash
# 1. 克隆
git clone https://github.com/songfy0118/diffprotacs-evaluator.git
cd diffprotacs-evaluator

# 2. 装依赖 (rdkit + numpy)
pip install -r requirements.txt

# 3. 跑 demo 验证装好了
python demo_data/make_demo_data.py
python -m evaluator.evaluate --xyz_dir demo_data/demo_output --out demo_eval.csv

# 4. 评估你自己的 DiffPROTACs 输出
python -m evaluator.evaluate --xyz_dir <你的_sample.py_输出目录> --out my_eval.csv
```

用 VS Code 打开整个目录,装一下 Python 扩展就能调试 / 看 CSV。
详细文档继续往下看。

---

## 快速开始 (3 步)

### 1. 安装依赖

```bash
pip install rdkit numpy
```

仅此而已。其他都是 Python 标准库。

### 2. 跑 demo (验证装好了)

```bash
cd diffprotacs_evaluator/
python demo_data/make_demo_data.py        # 生成 8 个分子 × 2 个采样的假 xyz
python -m evaluator.evaluate \
    --xyz_dir demo_data/demo_output \
    --out demo_eval.csv
```

应该看到类似这样的输出:
```
[evaluate] found 16 sampled molecules
[evaluate] validity rate (xyz -> mol): 14/16 = 87.50%

============================================================
DiffPROTACs Evaluation Summary
============================================================
  Total samples         : 16
  Valid (parses as mol) : 14  (87.50%)
  Unique SMILES         : 7  (50.00% of valid)
  Mean qed              : 0.5455
  Mean total_score      : 0.6505
  ...
```

并生成 `demo_eval.csv`,每个分子一行。

### 3. 评估你自己的 DiffPROTACs 输出

```bash
# 假设你跑过:
#   cd DiffPROTACs-main
#   python sample.py --xyz_path my_input.xyz --linker_len 8 --name target1 \
#                    --output_dir my_results/
# 然后:

python -m evaluator.evaluate \
    --xyz_dir /path/to/DiffPROTACs-main/my_results \
    --out my_eval.csv \
    --reference /path/to/protac_db_v3_smiles.txt    # 可选,用于 novelty
```

---

## 命令行参数

```
--xyz_dir          (必填) DiffPROTACs sample.py 的输出目录
                            (含 <uuid>/<i>_.xyz 子结构)
--out              输出 CSV 路径,默认 ./diffprotacs_eval.csv
--reference        参考 SMILES 文件 (.smi 或 .txt,一行一个),用于 novelty 打分
                   推荐: 训练集 / PROTAC-DB 3.0 / 已知专利分子库

# Docking 相关 (可选,不填则跳过 docking)
--docking_config   DockStream JSON 配置文件 (参考 configs/dockstream_template.json)
--dockstream_root  你的 DockStream-master-master 目录路径
--dockstream_python  跑 DockStream 的 Python 解释器,默认 "python"
                     (DockStream 依赖较重,通常需要单独的 conda env)
--dockstream_timeout 单次 docking 子进程超时秒数,默认 1800
```

---

## 接 DockStream / Vina 做真实 docking

1. 把 Protac-invent 那个 DockStream 装好(它需要 RDkit、AutoDockVina、相应受体 PDBQT 文件)
2. 复制 `configs/dockstream_template.json` 一份,把 5 个 `TODO_FILL_IN` 填掉:
   - Vina 二进制路径
   - 受体 PDBQT 路径
   - docking box 中心 (x,y,z) 和大小
   - 输出 poses SDF 路径
   - 输出 scores CSV 路径
3. 运行:

```bash
python -m evaluator.evaluate \
    --xyz_dir my_results \
    --out my_eval_with_docking.csv \
    --docking_config configs/my_target.json \
    --dockstream_root /path/to/DockStream-master-master
```

CSV 会多出 `raw_dock` (Vina 原始分,负值,越小越好) 和 `dock_score` ([0,1] 归一化,越大越好)两列。`total_score` 自动把 dock_score 计入。

> **小贴士**:Vina 打分慢,推荐先用纯 RDKit 模式快速过滤一遍,只对 `total_score` Top-N 的分子再上 docking。

---

## Python API 用法

详见 `examples/example_usage.py`。最常用三个动作:

```python
from evaluator import xyz_to_smiles, default_scorer, CompositeScorer
from evaluator import QEDScorer, SAScorer, SubstructureAlerts, NoveltyScorer

# 1. 单个 xyz 转 SMILES
smi = xyz_to_smiles("path/to/some.xyz")

# 2. 用默认打分器
scorer = default_scorer()
print(scorer.score_one(smi))
# -> {'validity': 1.0, 'qed': 0.59, 'sa_proxy': 0.25, ..., 'total_score': 0.65}

# 3. 自定义打分器组合
my_scorer = CompositeScorer(
    components=[
        (QEDScorer(),           1.0),
        (SAScorer(),            0.5),
        (SubstructureAlerts(),  2.0),       # 加重权重的硬门控
        (NoveltyScorer(my_ref), 1.0),
    ],
    mode="product",                          # 或 "sum"
)
```

---

## 文件结构

```
diffprotacs_evaluator/
├── README.md                          ← 你正在读的文件
├── evaluator/
│   ├── __init__.py
│   ├── xyz_to_smiles.py              ← 3D 坐标 → SMILES
│   ├── scorers.py                    ← 7 种 RDKit 打分器 + Composite
│   ├── docking_bridge.py             ← DockStream CLI 桥接(可选)
│   └── evaluate.py                   ← 主入口
├── configs/
│   └── dockstream_template.json      ← Vina docking 配置模板,填 5 个空就能用
├── demo_data/
│   ├── make_demo_data.py             ← 生成测试用 xyz(无需真正跑扩散模型)
│   └── reference_smiles.txt          ← novelty 测试用的小参考集
└── examples/
    └── example_usage.py              ← Python API 用法示例
```

---

## 一些技术说明 (重要,值得读)

### 1. xyz → SMILES 这一步是有损的

DiffPROTACs (以及所有 EDM 类扩散模型) 只输出**重原子的 3D 坐标**,不输出键级、不输出氢、不输出形式电荷。从 xyz 反推回合法分子是个老大难,文献中 EDM/MiDi/DiffSBDD 报告的 validity 普遍在 **30-90%** 之间,很大一部分波动就来自这一步。

本工具用的策略 (`evaluator/xyz_to_smiles.py`):
1. **距离阈值建初始键级** —— 用大量 PDB 统计的元素对距离窗口判断单/双/三键
2. **价态修复** —— 出现价态冲突时,迭代降低最长多键的键级
3. **芳香环修复** —— 检测被错判为累积双键的 6 元环,改写为 Kekulé 形式让 RDKit 重新芳香化
4. **rdDetermineBonds 回退** —— 几何路径失败时尝试 RDKit 官方 xyz2mol

**真实测试**(我们的 demo,8 类分子):

| 分子 | 转换结果 |
|---|---|
| benzene | ✅ `c1ccccc1` |
| toluene | ✅ `Cc1ccccc1` |
| naphthalene | ✅ `c1ccc2ccccc2c1` |
| anisole | ✅ `COc1ccccc1` |
| amide_linker (典型 PROTAC linker) | ✅ `CC(=O)NCCCCCNC(C)=O` |
| PEG 链 | ✅ `COCCOCCOCCOC` |
| caffeine (融合杂环) | ❌ (5+6 杂环组合对 xyz2mol 太难) |

> **更稳健的做法**(下一步可以做):用 DiffPROTACs 输入的 `frag_.xyz` 作为已知"锚",只对中间生成的 linker 部分做 xyz2mol。`evaluator/xyz_to_smiles.py` 里已经预留了 `reconstruct_with_fragment()` 函数,后续可以接进主流程。

### 2. 这套打分体系借鉴自 Protac-invent / REINVENT

Protac-invent 用 `reinvent_scoring` 的 `custom_product` 模式 (加权几何平均) 把多个打分组件融合。本工具的 `CompositeScorer(mode="product")` 完全等价,且**默认权重已经按 Protac-invent 的 BRD9 配置调好**。

主要差异:
- 我们用**纯 RDKit 重写**了打分组件,不依赖 `reinvent_scoring-0.0.73_bq.whl`,避免和 DiffPROTACs 的 conda env 冲突
- SA score 用了一个轻量代理而不是 RDKit Contrib 的 Ertl SA (避免额外依赖)。如果需要严格 Ertl SA,把 `sascorer.py` 从 RDKit Contrib 拷过来,改 `SAScorer.score()` 即可

### 3. 为什么 PROTAC 的 Lipinski 分数是软信号

PROTAC 分子设计上就违反 Ro5:
- 分子量通常 700-1200 (Ro5 ≤ 500)
- LogP 通常较高
- 氢键供受体多

所以本工具的 `LipinskiScorer` 返回的是"通过几条规则"的比例 (0~1),而不是硬切。在 PROTAC 评估里,**size_window** 和 **alerts_pass** 比 Lipinski 更值得关注。

### 4. 总分公式

`mode="product"` (默认):

$$\text{total} = \left(\prod_i s_i^{w_i}\right)^{1/\sum w_i}$$

—— 加权几何平均。任何一项是 0 (例如 invalid SMILES, 或命中 alert),总分立刻是 0。等价于 REINVENT 的 `custom_product`。

`mode="sum"`:

$$\text{total} = \frac{\sum_i w_i s_i}{\sum_i w_i}$$

—— 加权算术平均。任何一项是 0 不会拖死总分,适合"软"评估。

---

## 进一步可以做的(给 future-FF)

1. **anchor-based 评估**:利用 DiffPROTACs 的 `frag_.xyz` 作为已知部分,只对生成的 linker 做 SMILES 重建。骨架函数 `reconstruct_with_fragment()` 已经写好,接进 `evaluate.py` 主循环即可。
2. **接到 DiffPROTACs 训练循环里做 reward-guided sampling**:把 `CompositeScorer.score_one()` 当 reward function,在每个 diffusion step 用 reward gradient 引导采样。可以参考 DDPO (Black et al. 2023) 的实现。
3. **加 Ertl SA score**:从 RDKit Contrib 拷 `sascorer.py` 进 `evaluator/sa_score/`,替换 `SAScorer` 内部。
4. **加 3D 几何合理性检查**:键长偏离、原子相互穿透 (clash detection)。EDM 论文的 stability 指标。

---

## 致谢与参考

- **DiffPROTACs**: Bai Lab,基于 EDM-style equivariant diffusion
- **Protac-invent**: 基于 REINVENT 的 PROTAC RL 生成
- **DockStream**: MolecularAI/DockStream,docking 统一接口
- xyz2mol 算法: Kim & Kim, *J. Cheminform.* 2015
- QED: Bickerton et al., *Nat. Chem.* 2012

---

如有问题或需要扩展,记下来下次问 Claude。
