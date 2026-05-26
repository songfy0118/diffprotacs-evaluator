# 如何把这个评估器集成到 DiffPROTACs

## 选项 A: 作为独立工具用 (推荐,最简单)

直接把整个 `diffprotacs_evaluator/` 目录放在 `DiffPROTACs-main/` 旁边:

```
你的工作目录/
├── DiffPROTACs-main/         ← 朋友的原工程,不动
└── diffprotacs_evaluator/    ← 本工具
```

工作流:
```bash
# 1. 用 DiffPROTACs 生成
cd DiffPROTACs-main/
python sample.py \
    --xyz_path some_input.xyz \
    --linker_len 8 \
    --name brd9_test \
    --output_dir ../sample_outputs/run1

# 2. 切到 evaluator 目录评估
cd ../diffprotacs_evaluator/
python -m evaluator.evaluate \
    --xyz_dir ../sample_outputs/run1 \
    --out ../sample_outputs/run1_eval.csv

# 看结果
cat ../sample_outputs/run1_eval.csv | head
```

**好处**:DiffPROTACs 的 conda env 不需要装 evaluator 的依赖,反之亦然。两套环境完全独立,不会打架。

---

## 选项 B: 写一个 wrapper 一键跑

在 DiffPROTACs-main/ 根下放一个 `run_and_evaluate.sh`:

```bash
#!/bin/bash
# Usage: ./run_and_evaluate.sh <xyz_input> <linker_len> <name>

set -e

XYZ_INPUT=$1
LINKER_LEN=$2
NAME=$3
OUT_DIR="results_${NAME}_$(date +%Y%m%d_%H%M%S)"

# Activate DiffPROTACs env, sample
conda activate DiffPROTACs
python sample.py \
    --xyz_path "$XYZ_INPUT" \
    --linker_len "$LINKER_LEN" \
    --name "$NAME" \
    --output_dir "$OUT_DIR"

# Switch env, evaluate
conda activate diffprotacs_eval     # 一个轻量 env 只装了 rdkit + numpy
cd ../diffprotacs_evaluator
python -m evaluator.evaluate \
    --xyz_dir "../DiffPROTACs-main/$OUT_DIR" \
    --out "../DiffPROTACs-main/$OUT_DIR/eval.csv"

echo "Done. See ../DiffPROTACs-main/$OUT_DIR/eval.csv"
```

---

## 选项 C: 在 trainer.pred 内部直接打分 (进阶)

如果想让评估每生成一个 batch 就跑一次(例如训练时做 validation):

在 `DiffPROTACs-main/trainer.py` 的 `pred()` 方法最后加几行:

```python
# 在 trainer.py 末尾导入
import sys
sys.path.insert(0, "/path/to/diffprotacs_evaluator")
from evaluator import xyz_to_smiles, default_scorer

# 在 pred() 函数最后 (生成完所有 sample 后):
def pred(self, dataloader, output_dir, sample_fn=None):
    # ... 原有代码 ...

    # 新增:评估刚生成的所有 sample
    from evaluator.xyz_to_smiles import collect_xyz_files
    scorer = default_scorer()
    paths = collect_xyz_files(output_dir)
    n_valid, qed_sum = 0, 0.0
    for p in paths:
        smi = xyz_to_smiles(p)
        if smi:
            n_valid += 1
            qed_sum += scorer.score_one(smi)["qed"]
    print(f"[eval] validity: {n_valid}/{len(paths)}  "
          f"mean QED: {qed_sum / max(1, n_valid):.3f}")
```

> ⚠️ 这种用法 evaluator 必须装在 DiffPROTACs 同一个 conda env 里:
> ```bash
> conda activate DiffPROTACs
> pip install rdkit numpy   # DiffPROTACs 应该本来就装了 rdkit
> ```

---

## 选项 D: 把 evaluator 接到训练做 reward-guided sampling (最终目标)

这是朋友最终想做的:让评估器的分数反过来引导扩散过程,生成更高分的分子。

实现思路 (需要写代码,不是直接调用就行):

1. 在 `edm.py` 的去噪循环里,每若干步停下来用当前 `x_t` 解码出粗糙分子,跑 evaluator 算 reward
2. 用 reward 的梯度 (例如对 `x_t` 的 finite difference) 偏置下一步的去噪方向
3. 参考实现: **DDPO** (Black et al., *Diffusion Policy Optimization*, 2023) 或更简单的 **classifier-free guidance with reward**

这是一项**真正的研究工作**,不是几小时能做完的。给朋友的建议:
- 先用本 evaluator 做**离线评估**,跑出训练好的 DiffPROTACs 在不同条件下的 validity / QED / docking 分布
- 然后再考虑加 reward guidance —— 有了 baseline 才能证明 guidance 真的有提升

---

## 朋友最常会问的问题

**Q: validity 才 60% 是不是模型坏了?**
A: 这是 EDM 类扩散模型的**正常水平**。EDM 论文报告 ZINC 上 ~82%,QM9 上 ~92%,DrugMol 等更大数据上更低。如果你的 DiffPROTACs 在 PROTAC 数据上达到 50-70% validity,这就是合理的 baseline。

**Q: total_score 怎么看?**
A: 0 - 1 之间,越大越好。**>0.6** 是不错的分子,**>0.8** 是优质的(在默认 7 个组件的几何平均下,要每个组件都不差才能达到)。

**Q: 怎么提高 validity?**
A: 不是 evaluator 的问题,是模型本身的问题。常见方法:
1. 增大训练数据 / 训练步数
2. 用更严的几何约束 (例如 MiDi 的 bond-prediction head)
3. 后处理:对生成的 xyz 做短暂 MMFF 优化再评估 (键长收敛到化学合理范围)

**Q: 不想跑 docking 太慢,怎么办?**
A: 默认不跑 docking,只跑 RDKit 打分,3000 个分子大概 10 秒搞定。Docking 留到最后筛选 Top-100 的时候再上。

---

需要更详细的接入方案,问 Claude。
