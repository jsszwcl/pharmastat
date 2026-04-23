# PharmaStat - Statistical Visualization Software

Pharmaceutical data visualization and analysis tool with English UI and Chinese tooltip hints.

## Features

- **Group Comparison**: Bar+Error bars, Boxplot with automatic t-test significance markers
- **Dose-Response**: Sigmoid fitting, IC50/EC50 calculation
- **Correlation**: Scatter plot with regression line, r-value and p-value
- **Heatmap**: Data matrix heatmap visualization
- **Volcano Plot**: Differential expression analysis for genomics
- **Global Scroll**: Mouse wheel scrolling works across all tabs and input fields

## Visualization

- Matplotlib - Base plotting
- Seaborn - Statistical visualization with beautiful themes
- SciencePlots - Academic journal style

## Installation

```bash
pip install matplotlib pandas numpy scipy seaborn SciencePlots
```

## Usage

```bash
python pharmastat.py
```

## Dependencies

- Python 3.7+
- matplotlib
- pandas
- numpy
- scipy
- seaborn
- SciencePlots

## Interface

- All buttons and labels in English
- Hover over any element to see Chinese tooltip hints
- Group names: Control, Model, Low Dose, Mid Dose, High Dose, Positive