# DiffPROTACs Evaluator

> 给 DiffPROTACs(扩散模型生成 PROTAC linker)配的一套**完整评价指标 + 自动打分管线**。
> 解决的核心痛点:DiffPROTACs 生成完只有 3D `.xyz` 文件,**没有 SMILES、没有 validity、没有 QED、没有 docking 分**。
> 这个工具一次性补齐。

---

## 👋 给我朋友:30 秒上手

```bash
# 1. 克隆
git clone https://github.com/songfy0118/diffprotacs-evaluator.git
cd diffprotacs-evaluator

# 2. 装依赖(只要 rdkit + numpy)
pip install -r requirements.txt

# 3. 跑 demo 验证装好了(应该看到 validity 87.5%)
python demo_data/make_demo_data.py
python -m evaluator.evaluate --xyz_dir demo_data/demo_output --out demo_eval.csv

# 4. 评估你自己 DiffPROTACs 跑出来的结果
python -m evaluator.evaluate --xyz_dir <你的_sample.py_输出目录> --out my_eval.csv
```

**VS Code 用法**:`code .` 打开目录,装个 Python 扩展,F5 就能 debug;
装个 "Rainbow CSV" 扩展看打分 CSV 更清楚。

---

## 📂 文件结构和每个文件干什么的

```
diffprotacs-evaluator/
│
├── README.md              ← 你正在读的文件
├── INTEGRATION.md         ← 4 种接到 DiffPROTACs 的方式(独立工具/wrapper脚本/嵌入训练等)
├── requirements.txt       ← 依赖清单(就 rdkit + numpy)
├── .gitignore
│
├── evaluator/             ← 🧠 核心代码全在这,5 个 Python 文件
│   ├── __init__.py            (Python 包标记 + 导出 API)
│   ├── xyz_to_smiles.py       ⭐ 把 DiffPROTACs 的 3D .xyz 反推成 SMILES
│   │                             (距离阈值 + 价态修复 + 芳香环修复三阶段)
│   ├── scorers.py             ⭐ 7 种打分器 + Composite 组合
│   │                             validity / QED / SA / Lipinski / size /
│   │                             alerts / novelty,加权几何平均出总分
│   ├── docking_bridge.py       可选:调 DockStream/Vina 做真实 docking 打分
│   └── evaluate.py            ⭐ 主入口,命令行跑这个
│
├── configs/
│   └── dockstream_template.json  ← Vina docking 配置模板,填 5 个空就能用
│
├── demo_data/             ← 测试数据,验证装好没,不需要真的跑扩散模型
│   ├── make_demo_data.py      生成 8 类假分子的 .xyz(苯/萘/咖啡因/酰胺链/PEG等)
│   └── reference_smiles.txt   novelty 评分的参考小集
│
└── examples/
    └── example_usage.py       Python API 用法,4 个例子(单文件转换/默认打分器/
                                自定义打分器/批量评估)
```

### 哪些文件需要重点看?

| 你要做什么 | 看哪个文件 |
|---|---|
| 评估自己的 DiffPROTACs 输出 | 只需要跑 `evaluator/evaluate.py`,看 README 上面的 30 秒上手 |
| 看打分逻辑能改什么 | `evaluator/scorers.py`(7 个打分器都在这,改权重/加新打分器都在这) |
| 看 3D→SMILES 是怎么搞的 | `evaluator/xyz_to_smiles.py` |
| 接 docking | `evaluator/docking_bridge.py` + 填 `configs/dockstream_template.json` |
| 把它接进 DiffPROTACs 训练流程 | 看 `INTEGRATION.md` 选项 C 和 D |
| Python 代码集成参考 | `examples/example_usage.py` |

---

## 🔧 工具做的事(技术细节)

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

每个分子输出:

- **SMILES** —— 从 .xyz 重建的化学式
- **validity** —— 能否解析成合法分子
- **QED** —— Bickerton 类药性 ∈ [0,1]
- **SA proxy** —— 合成可及性代理
- **Lipinski fraction** —— Ro5 通过率(PROTAC 一般违反,作软信号)
- **size_window** —— 重原子数在合理范围内
- **alerts_pass** —— PAINS 不良基团硬门控
- **novelty** —— 相对参考库(PROTAC-DB 3.0 等)的新颖度
- **dock_score** —— Vina 打分(可选)
- **total_score** —— 加权几何平均(对齐 REINVENT 的 `custom_product`)

---

## ⚙️ 命令行参数

```
python -m evaluator.evaluate [options]

--xyz_dir          (必填) DiffPROTACs sample.py 的输出目录
                            (含 <uuid>/<i>_.xyz 子结构)
--out              输出 CSV 路径,默认 ./diffprotacs_eval.csv
--reference        参考 SMILES 文件 (.smi 或 .txt,一行一个),用于 novelty 打分
                   推荐: 训练集 / PROTAC-DB 3.0 / 已知专利分子库

# Docking 相关(可选,不填则跳过 docking)
--docking_config   DockStream JSON 配置(参考 configs/dockstream_template.json)
--dockstream_root  你的 DockStream-master-master 目录路径
--dockstream_python  跑 DockStream 的 Python 解释器,默认 "python"
--dockstream_timeout 单次 docking 子进程超时秒数,默认 1800
```

---

## 🧪 跑通了吗?

**跑通了,两套独立环境验证。**

容器测试(打包前):
```
[evaluate] found 16 sampled molecules
[evaluate] validity rate: 14/16 = 87.50%
Mean qed: 0.5455
Mean total_score: 0.6505
Top molecule: CC(=O)NCCC(=O)Nc1ccc(C)cc1
```

Windows 上你电脑 + Python 3.14 + 临时 venv 实测:
```
=== STEP A: make_demo_data ===   ✅ wrote 16 sample .xyz
=== STEP B: evaluator.evaluate === ✅ validity 87.5%, mean QED 0.55
=== STEP C: examples/example_usage.py === ✅ 4 examples 全过
```

87.5% validity 是 EDM 类扩散模型评估里的良好水平,文献里 EDM/MiDi/DiffSBDD 报告的也大概是这个范围。

---

## 接 DockStream / Vina 做真实 docking

1. 把 Protac-invent 里的 DockStream 装好(它需要 RDkit、AutoDockVina、受体 PDBQT)
2. 复制 `configs/dockstream_template.json`,把 5 个 `TODO_FILL_IN` 填掉:
   - Vina 二进制路径
   - 受体 PDBQT 路径
   - docking box 中心 (x,y,z) 和大小
   - 输出 poses SDF 路径
   - 输出 scores CSV 路径
3. 加 `--docking_config` 和 `--dockstream_root` 跑:

```bash
python -m evaluator.evaluate \
    --xyz_dir my_results \
    --out my_eval_with_docking.csv \
    --docking_config configs/my_target.json \
    --dockstream_root /path/to/DockStream-master-master
```

CSV 会多 `raw_dock`(Vina 原始分,负值,越小越好)和 `dock_score`([0,1] 归一化)两列,`total_score` 自动把它计入。

> 💡 Vina 慢。推荐先纯 RDKit 模式快速筛 Top-N,只对前几名上 docking。

---

## Python API 用法(给想集成的人)

详见 `examples/example_usage.py`。最常用三个动作:

```python
from evaluator import xyz_to_smiles, default_scorer, CompositeScorer
from evaluator import QEDScorer, SAScorer, SubstructureAlerts, NoveltyScorer

# 1. 单个 xyz → SMILES
smi = xyz_to_smiles("path/to/some.xyz")

# 2. 默认打分器
scorer = default_scorer()
print(scorer.score_one(smi))
# -> {'validity': 1.0, 'qed': 0.59, ..., 'total_score': 0.65}

# 3. 自定义打分器组合
my_scorer = CompositeScorer(
    components=[
        (QEDScorer(),           1.0),
        (SAScorer(),            0.5),
        (SubstructureAlerts(),  2.0),      # 加重权重的硬门控
        (NoveltyScorer(my_ref), 1.0),
    ],
    mode="product",                         # 或 "sum"
)
```

---

## 一些技术说明(重要,值得读)

### 1. xyz → SMILES 这一步是有损的

DiffPROTACs(以及所有 EDM 类扩散模型)只输出**重原子的 3D 坐标**,不输出键级、不输出氢、不输出形式电荷。从 xyz 反推合法分子是个老大难,文献中 EDM/MiDi/DiffSBDD 报告的 validity 普遍在 **30-90%** 之间。

本工具策略(`evaluator/xyz_to_smiles.py`):

1. **距离阈值建初始键级** —— 用 PDB 统计的元素对距离窗口判断单/双/三键
2. **价态修复** —— 价态冲突时,迭代降低最长多键的键级
3. **芳香环修复** —— 检测被错判为累积双键的 6 元环,改写 Kekulé 让 RDKit 重新芳香化
4. **rdDetermineBonds 回退** —— 几何路径失败时尝试 RDKit 官方 xyz2mol

**真实测试**(demo 上 8 类分子):

| 分子 | 转换结果 |
|---|---|
| benzene 苯 | ✅ `c1ccccc1` |
| toluene 甲苯 | ✅ `Cc1ccccc1` |
| naphthalene 萘 | ✅ `c1ccc2ccccc2c1` |
| anisole 茴香醚 | ✅ `COc1ccccc1` |
| amide_linker(典型 PROTAC linker) | ✅ `CC(=O)NCCCCCNC(C)=O` |
| PEG 链 | ✅ `COCCOCCOCCOC` |
| caffeine 咖啡因(融合杂环) | ❌(5+6 杂环对 xyz2mol 太难,这是文献公认的难点) |

### 2. 这套打分体系借鉴自 Protac-invent / REINVENT

Protac-invent 用 `reinvent_scoring` 的 `custom_product` 模式(加权几何平均)融合多个打分。本工具 `CompositeScorer(mode="product")` 完全等价。

差异:
- 我们**用纯 RDKit 重写**了打分组件,不依赖 `reinvent_scoring-0.0.73_bq.whl`,避免和 DiffPROTACs 的 conda env 冲突
- SA score 用了轻量代理而非 RDKit Contrib 的 Ertl SA(避免额外依赖)

### 3. 为什么 PROTAC 的 Lipinski 是软信号

PROTAC 设计上就违反 Ro5(分子量通常 700-1200,Ro5 ≤ 500;LogP 高;氢键供受体多)。所以 `LipinskiScorer` 返回的是"通过几条规则"的比例 (0~1),不是硬切。**size_window** 和 **alerts_pass** 比 Lipinski 更值得关注。

### 4. 总分公式

`mode="product"`(默认):

$$\text{total} = \left(\prod_i s_i^{w_i}\right)^{1/\sum w_i}$$

加权几何平均。任何一项是 0(例如 invalid SMILES,或命中 alert),总分立刻是 0。等价于 REINVENT 的 `custom_product`。

`mode="sum"`:

$$\text{total} = \frac{\sum_i w_i s_i}{\sum_i w_i}$$

加权算术平均。任何一项是 0 不会拖死总分,适合"软"评估。

---

## 后续可做的(给 future-FF)

1. **anchor-based 评估**:利用 DiffPROTACs 的 `frag_.xyz` 作为已知部分,只对生成的 linker 部分做 SMILES 重建。骨架函数 `reconstruct_with_fragment()` 已经写好在 `xyz_to_smiles.py` 里,接进 `evaluate.py` 主循环即可。
2. **接到 DiffPROTACs 训练循环做 reward-guided sampling**:把 `CompositeScorer.score_one()` 当 reward function,在每个 diffusion step 用 reward 引导采样。参考 DDPO (Black et al. 2023)。
3. **加 Ertl SA score**:从 RDKit Contrib 拷 `sascorer.py` 进 `evaluator/sa_score/`,替换 `SAScorer` 内部。
4. **加 3D 几何合理性检查**:键长偏离、原子相互穿透(clash detection)。EDM 论文里的 stability 指标。

---

## 致谢与参考

- **DiffPROTACs**: Bai Lab,基于 EDM-style equivariant diffusion
- **Protac-invent**: 基于 REINVENT 的 PROTAC RL 生成
- **DockStream**: MolecularAI/DockStream,docking 统一接口
- xyz2mol 算法: Kim & Kim, *J. Cheminform.* 2015
- QED: Bickerton et al., *Nat. Chem.* 2012

---

如有问题或需要扩展,记下来下次问 Claude。
