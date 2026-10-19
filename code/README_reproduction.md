# 统计复现入口与本轮核验范围

本文件给出当前包内历史统计入口的实际位置、必要参数和已记录的依赖版本。2026-10-17 本轮只做输入/路径/指纹静态核验、既有结果比对及稿件图件整理；**未重拟合回归、未运行 wild bootstrap**。下面的拟合命令仅供以后在备妥环境时复算使用，不是本轮执行记录。

## 先做无拟合预检

在本包根目录运行：

```sh
python3 code/check_reproduction_inputs.py
```

预检只用 Python 标准库，校验问卷和历史代码/参考结果的 SHA-256、入口实际解析的包根目录、追加诊断所需的初始结果路径，以及活动学校缩放表的范围。它不会导入统计分析脚本、拟合模型、枚举 bootstrap 或安装软件。默认仅将 JSON 写到终端；要保留报告，可增加 `--report /绝对路径/已存在目录/新的预检报告.json`，该文件必须位于包外且尚不存在。

预检分别报告输入是否通过、当前 Python 能否找到所需软件包，以及历史 audit 记录的版本。输入通过不代表数值复现成功，也不代表当前 Python 已装齐依赖。本轮指定的 bundle Python 包元数据中，所列统计依赖仅检出 NumPy/pandas；Patsy、SciPy、statsmodels、scikit-learn、threadpoolctl 未被检出，因而本轮没有进行统计执行验证。没有尝试安装或以其他环境绕过此限制。

## 可定位的历史入口

只把当前包第一层 `archive/原临床分析代码/` 下的脚本作为以下命令的入口。上一版完整快照若再次放入 archive，其额外目录深度会改变历史 `Path.parents` 的解析；该快照用于存档，不用作可运行入口。

| 历史任务 | 本包入口 | 主要输入/参考结果 | 已有输出参照 |
|---|---|---|---|
| 12 October 内容 common-prior、四子量表与学段描述 | `archive/原临床分析代码/revision_20261012/content_analysis.py` | `data/questionnaire_analysis.csv`、同目录 `reference_results/`、`CONTENT_PLAN.md`、上级 `revision_common.py`/`revision_runtime.py` | `archive/原分析汇总结果/revision_20261012/` |
| 12 October 留存/IPW诊断 | `archive/原临床分析代码/revision_20261012/ipw_analysis.py` | 同一问卷、同目录参考结果、`IPW_PLAN.md`、上级输出保护代码 | 同上 |
| 15 October 年级 common-prior、函数形式/影响、协方差与EN−PN | `archive/原临床分析代码/post_review_20261015/post_review_analysis.py` | 同一问卷、`dependencies/content_analysis_20261012.py`、`reference_results/`、`POST_REVIEW_PLAN_20261015.md`、上级共同代码 | `archive/原分析汇总结果/revision_20261015/` |
| 15 October 两个混合规格、既有全样条重建与曲线 | `archive/原临床分析代码/post_review_20261015/diagnostic_addendum.py` | 上一行的代码与参考输入；显式 `--initial-results` 指向已有初始结果；`POST_REVIEW_DIAGNOSTIC_ADDENDUM_20261015.md` | `archive/原分析汇总结果/diagnostic_addendum_20261015/` |

上述4个入口的 `Path(__file__).resolve().parents[3]` 在当前布局解析为本包根目录。post-review 脚本加载的依赖由其入口设置正确的问卷与参考位置；不要把 `dependencies/content_analysis_20261012.py` 当独立命令运行。

历史源文件保留不改。所有拟合入口要求 `--output` 为**包外、此前不存在的新目录**，以保留冻结结果。工作目录设为本包根目录后，以下是完整参数示例；把 `/绝对路径/新的复现输出/...` 替换为计划使用的新包外位置。各命令互相独立；追加诊断可直接使用已保留的初始结果，无需为它重新运行前一分析。

```sh
python3 archive/原临床分析代码/revision_20261012/content_analysis.py \
  --output /绝对路径/新的复现输出/content_20261012

python3 archive/原临床分析代码/revision_20261012/ipw_analysis.py \
  --output /绝对路径/新的复现输出/ipw_20261012

python3 archive/原临床分析代码/post_review_20261015/post_review_analysis.py \
  --output /绝对路径/新的复现输出/post_review_20261015

python3 archive/原临床分析代码/post_review_20261015/diagnostic_addendum.py \
  --initial-results archive/原分析汇总结果/revision_20261015 \
  --output /绝对路径/新的复现输出/diagnostic_addendum_20261015
```

`diagnostic_addendum.py` 保留的默认 `paper/evidence/results/revision_20261015` 是原包布局，在当前包中不存在，因此上例的 `--initial-results` **必须显式提供**。目标目录内已有 `spline_knots.json` 与 `function_form_iqr_contrasts.csv`。若另行使用归档 `build_diagnostic_figure.py`，同样要显式给 `--source archive/原分析汇总结果/diagnostic_addendum_20261015` 与新的包外 `--output`；投稿活动图件应使用当前 `code/build_current_figures.py`，详见 `README_figures.md`。

输出CSV如需比较，可使用 `archive/原临床分析代码/compare_revision_outputs.py --actual <新的输出目录> --frozen <相应既有结果目录> --report <新的JSON路径>`；它比较实际目录中的聚合CSV，不执行模型。输出集合须与所选冻结目录相匹配，不能把不同阶段或不同规格混为同一复算。

## 历史拟合环境与本轮整理环境

下列数字直接来自历史统计 audit JSON，不以本轮绘图或 Word 环境版本代替：

| 依赖 | 历史记录版本 | 记录来源 |
|---|---|---|
| Python | 3.12.14 | 12 October 内容/IPW与15 October首轮 audit |
| NumPy | 2.3.5 | 上述 audit及追加诊断 audit |
| pandas | 2.2.3 | 同上 |
| SciPy | 1.18.1 | 同上 |
| statsmodels | 0.15.0 | 同上 |
| Patsy | 1.0.3 | 同上 |
| scikit-learn（导入名sklearn） | 1.9.1 | `archive/原分析汇总结果/revision_20261012/ipw_summary.json` |
| threadpoolctl | 历史 audit 未记录版本 | 历史脚本直接导入；复算环境必须提供，不能虚构其历史版本 |

版本来源文件为：

- `archive/原分析汇总结果/revision_20261012/content_revision_audit.json`
- `archive/原分析汇总结果/revision_20261012/ipw_summary.json`
- `archive/原分析汇总结果/revision_20261015/post_review_analysis_audit.json`
- `archive/原分析汇总结果/diagnostic_addendum_20261015/diagnostic_addendum_audit.json`

最后一份追加诊断 audit 记录了 NumPy/pandas/SciPy/statsmodels/Patsy版本，但没有Python版本字段；上表不把缺少的字段解释为已由该文件核实。以上为已记录版本，不是经过本轮安装验证的锁定环境文件。

`code/runtime_versions.json` 的 document/figure 字段记录本轮Word/打包环境和本轮图件实际导入的依赖。这些元数据与历史统计拟合环境各有范围，即使部分数字相同，也不能据此声称本轮已重跑统计模型。

## 结果与数据的范围

- `source_data/analysis/core_26_model_catalog.csv` 是12 October已完成的26项核心模型计算核验登记与来源索引；本轮保留历史核验记录，未重跑26项。以上4个历史入口不能被描述为本轮26项全套重拟合或一键复现证明。
- 活动 `source_data/analysis/scaling_register_full.csv` 仅含学校队列的24行。波次1–2使用完整合格W1基线N=3,342，波次2–3和其限制样本使用完整W2基线N=8,902；标准差ddof=0。
- 归档的历史15 October流程还会读取 `reference_results/external_clinical_scaling.json` 和 `external_clinical_sample.csv`，以保留其历史缩放登记输出；这些文件存在，但成人常数不属于本稿活动分析。若以后重跑，历史 `scaling_register_full.csv` 可能仍含成人行，不能直接覆盖活动学校版。成人原始个体数据未在本轮读取或核实。
- `questionnaire_analysis.csv` 是去标识的处理后问卷表，不是原始题本或私有原始导出。输入/计算核验不能审定CTQ原始计分、实际施测、伦理、同意或数据共享许可。
- 当前绘图入口沿用原投稿包的聚合结果/原图依赖，必须显式给 `--source-package`；本说明未把它改称无需原包的自包含绘图入口。具体命令和本轮图件修订由 `README_figures.md` 说明。
