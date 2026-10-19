# 图件代码

`build_current_figures.py` 只读取既有聚合CSV并准备活动图件及其源表，不执行模型拟合，不导入或运行历史模型入口。

运行需要 Python、NumPy、pandas、matplotlib，以及 `figure_qa_tools/audit_panel_alignment.py`。四个核查脚本从本机安装的 nature-figure 技能原样保留；PDF碰撞检查另需PyMuPDF。命令参数使用明确的输入和新输出目录：

```text
python build_current_figures.py --source-package <已解压的原投稿包> --output <新的输出目录> --qa <新的QA目录> --audit-scripts <核查工具目录>
```

`--source-package` 用于读取保存的估计和复制三张未改原图。实际活动panel CSV及其原始相对位置/哈希已列入 `source_data/figure_source_registry.csv`；不依赖私有问卷原始导出或录音。Figure2使用既有留存、组成与prior-distress区间估计；Figure3使用既有common-prior及描述资料；S3只修复Q75辅助线穿图例问题，不改变任何曲线、区间或频数。

当前三张原图1/S1/S2按字节保留；2/3/S3按PNG/PDF/SVG导出，PNG600dpi。所有新绘制图均通过rendered panel alignment及PDF碰撞核查。旧图1/S1的刻意图内注释叠淡背景没有获得自动collision PASS，见单独QA记录；不应声称六图全部自动通过。

2026-10-17 收尾修复沿用同一 Python/matplotlib 绘图流程。`fig1c_main_distress_units.csv` 和 `figS2_reference_main_estimates.csv` 的重建明确筛选 `sample=All`，各从原六行保留两行，并原样保留选中 CSV 行的数字字符串及空值；完整子集结果继续保存在 archive。

Figure 2 的真实条形斜线仍按原方法栅格化。图例斜线按整个 DrawingArea 栅格化，并把容器偏移同步到切换后的栅格 DPI 坐标，以免单独栅格化代理 Patch 将其错误导出到左下边距。图例文字继续是矢量文本。此次 Figure 2 PNG 与上一版逐字节一致，仅 PDF/SVG 的错误图例位置得到修复。

可加 `--only-figure2` 在独立 staging 目录准备全部源数据映射并仅渲染 Figure 2，其他五图不会写出。运行后按 scope 精确复制两份源表、registry 与 Figure 2 文件，不把 staging 中的其他数据覆盖到交付目录。

本轮实际导入依赖经过验证：指定 primary-runtime Python 提供 NumPy 2.3.5、pandas 2.2.3、Pillow 12.3.0；原分析环境 site-packages 以追加路径提供 matplotlib 3.11.2 与 PyMuPDF 1.28.2。实际路径及命令见 `validation/figure_cleanup_F01_F04/runtime_dependencies.json`；版本记录不代替导入检查。通过 fresh render alignment、PDF 6.0 pt 最小文字和 strict collision 门禁，源代码预检为20 PASS/0 FAIL，唯一 WARN 是当前授权导出没有 TIFF，保留既定 PNG/PDF/SVG 三格式。
