# PharmaStat - 药学可视化统计软件

药学数据可视化分析工具，支持柱状图、箱线图、量效曲线、相关性分析、热图等。

## 功能

- **组间比较**: 柱状图+误差线、箱线图，自动t检验显著性标注
- **量效曲线**: S型拟合，计算IC50/EC50
- **相关性**: 散点图+回归直线，计算r值和p值
- **热图**: 数据矩阵热力图可视化

## 安装

```bash
pip install matplotlib pandas numpy scipy
```

## 使用

```bash
python pharmastat.py
```

## 依赖

- Python 3.7+
- matplotlib
- pandas
- numpy
- scipy