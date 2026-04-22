#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PharmaStat - 药学可视化的统计软件
支持柱状图+误差线、箱线图、量效曲线、相关性热图等
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

import seaborn as sns
sns.set_theme(style="whitegrid", palette="muted", font_scale=1.1)

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

try:
    plt.style.use('science')
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'Arial Unicode MS', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
except:
    pass

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

GROUP_ORDER = ['Control', 'Model', 'Low Dose', 'Mid Dose', 'High Dose', 'Positive']


class VolcanoGrid(ttk.Frame):
    def __init__(self, parent, rows=20):
        super().__init__(parent, relief=tk.SUNKEN, borderwidth=1)
        self.rows = rows
        self.cells = {}
        self.setup_grid()
        
    def setup_grid(self):
        headers = ['Gene Name', 'Log2FC', 'P-value']
        for c, h in enumerate(headers):
            ttk.Label(self, text=h, width=15, anchor='center', font=('Arial', 10, 'bold')).grid(row=0, column=c, padx=2, pady=3)
            
        for r in range(self.rows):
            for c in range(3):
                entry = ttk.Entry(self, width=15, justify='center')
                entry.grid(row=r+1, column=c, padx=2, pady=2)
                self.cells[(r, c)] = entry
                
    def add_row(self):
        if self.rows >= 100:
            return
        self.rows += 1
        r = self.rows - 1
        for c in range(3):
            entry = ttk.Entry(self, width=15, justify='center')
            entry.grid(row=r+1, column=c, padx=2, pady=2)
            self.cells[(r, c)] = entry
            
    def del_row(self):
        if self.rows <= 1:
            return
        r = self.rows - 1
        for c in range(3):
            self.cells[(r, c)].grid_forget()
            self.cells[(r, c)].destroy()
            del self.cells[(r, c)]
        self.rows -= 1
        
    def get_data(self):
        genes = []
        log2fcs = []
        pvals = []
        for r in range(self.rows):
            gene = self.cells[(r, 0)].get().strip()
            log2fc = self.cells[(r, 1)].get().strip()
            pval = self.cells[(r, 2)].get().strip()
            
            if gene:
                try:
                    genes.append(gene)
                    log2fcs.append(float(log2fc) if log2fc else 0)
                    pvals.append(float(pval) if pval else 1)
                except ValueError:
                    pass
        return genes, log2fcs, pvals
GROUP_COLORS = ['#4CAF50', '#FF5722', '#2196F3', '#9C27B0', '#FF9800', '#795548']
DEFAULT_INDICATORS = ['血压(mmHg)', '血糖(mmol/L)', '转氨酶(U/L)', '脏器系数(%)', '痛阈值(s)', '收缩力(mN)']

GROUP_TOOLTIP = {
    'Control': '对照组',
    'Model': '模型组',
    'Low Dose': '低剂量组',
    'Mid Dose': '中剂量组',
    'High Dose': '高剂量组',
    'Positive': '阳性药组',
}

TOOLTIP_MAP = {
    'Group': '组间比较',
    'Dose-Response': '量效曲线',
    'Correlation': '相关性分析',
    'Heatmap': '热图',
    'Indicator': '指标',
    'Clear': '清空数据',
    '+ Add Col': '添加样本列',
    '- Del Col': '删除样本列',
    'Bar+Error': '柱状图+误差线',
    'Boxplot': '箱线图',
    'Import CSV': '导入CSV文件',
    'Export': '导出统计结果',
    '+ Add Row': '添加一行',
    '- Del Row': '删除一行',
    'Generate Curve': '生成量效曲线',
    'Scatter Plot': '散点图+回归',
    '+ Add Row': '添加行',
    '- Del Row': '删除行',
    '+ Add Col': '添加列',
    '- Del Col': '删除列',
    'Generate': '生成热图',
    'Save PNG': '保存为PNG图片',
    'Save PDF': '保存为PDF文档',
}

def add_tooltip(widget, text_key):
    tooltip_text = TOOLTIP_MAP.get(text_key, '')
    ToolTip(widget, tooltip_text)

class ToolTip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tipwindow = None
        self.widget.bind('<Enter>', self.show)
        self.widget.bind('<Leave>', self.hide)
        
    def show(self, event=None):
        if self.tipwindow or not self.text:
            return
        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 5
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_geometry(f'+{x}+{y}')
        tw.wm_overrideredirect(True)
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                      font=('Microsoft YaHei', 10), bg='#ffffe0', fg='#000000',
                      relief=tk.SOLID, borderwidth=1, padx=5, pady=3)
        label.pack()
        
    def hide(self, event=None):
        if self.tipwindow:
            self.tipwindow.destroy()
            self.tipwindow = None


class DataGrid(ttk.Frame):
    def __init__(self, parent, rows=6, cols=8):
        super().__init__(parent, relief=tk.SUNKEN, borderwidth=1)
        self.rows = rows
        self.cols = cols
        self.cells = {}
        self.setup_grid()
        
    def setup_grid(self):
        ttk.Label(self, text='Sample→', font=('Arial', 9), width=8).grid(row=0, column=0, padx=2, pady=3)
        
        for c in range(self.cols):
            label = ttk.Label(self, text=f"#{c+1}", width=6, anchor='center', font=('Arial', 9))
            label.grid(row=0, column=c+1, padx=1, pady=3)
            
        for r in range(self.rows):
            label_text = GROUP_ORDER[r]
            label = ttk.Label(self, text=label_text, width=10, anchor='w')
            label.grid(row=r+1, column=0, padx=2, pady=2)
            if label_text in GROUP_TOOLTIP:
                ToolTip(label, GROUP_TOOLTIP[label_text])
            for c in range(self.cols):
                entry = ttk.Entry(self, width=6, justify='center')
                entry.grid(row=r+1, column=c+1, padx=1, pady=2)
                self.cells[(r, c)] = entry
                
    def add_col(self):
        if self.cols >= 100:
            return
        self.cols += 1
        c = self.cols - 1
        label = ttk.Label(self, text=f"#{c+1}", width=6, anchor='center', font=('Arial', 9))
        label.grid(row=0, column=c+1, padx=1, pady=3)
        for r in range(self.rows):
            entry = ttk.Entry(self, width=6, justify='center')
            entry.grid(row=r+1, column=c+1, padx=1, pady=2)
            self.cells[(r, c)] = entry
            
    def del_col(self):
        if self.cols <= 1:
            return
        c = self.cols - 1
        for r in range(self.rows):
            self.cells[(r, c)].grid_forget()
            self.cells[(r, c)].destroy()
            del self.cells[(r, c)]
        self.cols -= 1
        
    def get_data(self):
        data = {}
        for r, group in enumerate(GROUP_ORDER):
            values = []
            for c in range(self.cols):
                val = self.cells[(r, c)].get().strip()
                if val:
                    try:
                        values.append(float(val))
                    except ValueError:
                        pass
            if values:
                data[group] = values
        return data
        
    def set_data(self, data):
        for r in range(self.rows):
            for c in range(self.cols):
                self.cells[(r, c)].delete(0, tk.END)
        for group, vals in data.items():
            if group in GROUP_ORDER:
                r = GROUP_ORDER.index(group)
                for c, val in enumerate(vals):
                    if c < self.cols:
                        self.cells[(r, c)].insert(0, str(val))


class DoseResponseGrid(ttk.Frame):
    def __init__(self, parent, rows=8):
        super().__init__(parent, relief=tk.SUNKEN, borderwidth=1)
        self.rows = rows
        self.cells = {}
        self.setup_grid()
        
    def setup_grid(self):
        headers = ['Dose', 'Effect Mean', 'Effect SD', 'N']
        for c, h in enumerate(headers):
            ttk.Label(self, text=h, width=12, anchor='center', font=('Arial', 10, 'bold')).grid(row=0, column=c, padx=2, pady=3)
            
        for r in range(self.rows):
            for c in range(4):
                entry = ttk.Entry(self, width=12, justify='center')
                entry.grid(row=r+1, column=c, padx=2, pady=2)
                self.cells[(r, c)] = entry
                
    def add_row(self):
        if self.rows >= 50:
            return
        self.rows += 1
        r = self.rows - 1
        for c in range(4):
            entry = ttk.Entry(self, width=12, justify='center')
            entry.grid(row=r+1, column=c, padx=2, pady=2)
            self.cells[(r, c)] = entry
            
    def del_row(self):
        if self.rows <= 1:
            return
        r = self.rows - 1
        for c in range(4):
            self.cells[(r, c)].grid_forget()
            self.cells[(r, c)].destroy()
            del self.cells[(r, c)]
        self.rows -= 1
        
    def get_data(self):
        doses = []
        effects = []
        sds = []
        ns = []
        for r in range(self.rows):
            dose = self.cells[(r, 0)].get().strip()
            effect = self.cells[(r, 1)].get().strip()
            sd = self.cells[(r, 2)].get().strip()
            n = self.cells[(r, 3)].get().strip()
            
            if dose and effect:
                try:
                    doses.append(float(dose))
                    effects.append(float(effect))
                    sds.append(float(sd) if sd else 0)
                    ns.append(int(n) if n else 1)
                except ValueError:
                    pass
        return doses, effects, sds, ns
        
    def set_data(self, doses, effects, sds, ns):
        for r in range(self.rows):
            for c in range(4):
                self.cells[(r, c)].delete(0, tk.END)
                
        for i, (d, e, s, n) in enumerate(zip(doses, effects, sds, ns)):
            self.cells[(i, 0)].insert(0, str(d))
            self.cells[(i, 1)].insert(0, str(e))
            self.cells[(i, 2)].insert(0, str(s))
            self.cells[(i, 3)].insert(0, str(n))


class CorrelationGrid(ttk.Frame):
    def __init__(self, parent, max_rows=30):
        super().__init__(parent, relief=tk.SUNKEN, borderwidth=1)
        self.max_rows = max_rows
        self.cells = {}
        self.setup_grid()
        
    def setup_grid(self):
        ttk.Label(self, text='X Variable', width=15, anchor='center', font=('Arial', 10, 'bold')).grid(row=0, column=0, padx=2, pady=3)
        ttk.Label(self, text='Y Variable', width=15, anchor='center', font=('Arial', 10, 'bold')).grid(row=0, column=1, padx=2, pady=3)
        ttk.Label(self, text='Group (Opt)', width=12, anchor='center', font=('Arial', 10, 'bold')).grid(row=0, column=2, padx=2, pady=3)
            
        for r in range(self.max_rows):
            for c in range(3):
                entry = ttk.Entry(self, width={0: 15, 1: 15, 2: 12}[c], justify='center')
                entry.grid(row=r+1, column=c, padx=2, pady=1)
                self.cells[(r, c)] = entry
                
    def add_row(self):
        if self.max_rows >= 200:
            return
        self.max_rows += 1
        r = self.max_rows - 1
        for c in range(3):
            entry = ttk.Entry(self, width={0: 15, 1: 15, 2: 12}[c], justify='center')
            entry.grid(row=r+1, column=c, padx=2, pady=1)
            self.cells[(r, c)] = entry
            
    def del_row(self):
        if self.max_rows <= 2:
            return
        r = self.max_rows - 1
        for c in range(3):
            self.cells[(r, c)].grid_forget()
            self.cells[(r, c)].destroy()
            del self.cells[(r, c)]
        self.max_rows -= 1
                
    def get_data(self):
        x_data = []
        y_data = []
        groups = []
        for r in range(self.max_rows):
            x = self.cells[(r, 0)].get().strip()
            y = self.cells[(r, 1)].get().strip()
            g = self.cells[(r, 2)].get().strip()
            
            if x and y:
                try:
                    x_data.append(float(x))
                    y_data.append(float(y))
                    groups.append(g if g else 'Default')
                except ValueError:
                    pass
        return x_data, y_data, groups
        
    def set_data(self, x_data, y_data, groups):
        for r in range(self.max_rows):
            for c in range(3):
                self.cells[(r, c)].delete(0, tk.END)
                
        for i, (x, y, g) in enumerate(zip(x_data, y_data, groups)):
            self.cells[(i, 0)].insert(0, str(x))
            self.cells[(i, 1)].insert(0, str(y))
            if g:
                self.cells[(i, 2)].insert(0, str(g))


class HeatmapGrid(ttk.Frame):
    def __init__(self, parent, rows=10, cols=10):
        super().__init__(parent, relief=tk.SUNKEN, borderwidth=1)
        self.rows = rows
        self.cols = cols
        self.cells = {}
        self.setup_grid()
        
    def setup_grid(self):
        ttk.Label(self, text='  Heatmap Data Matrix  ', font=('Arial', 10, 'bold')).grid(row=0, column=0, columnspan=self.cols+1, pady=5)
            
        for r in range(self.rows):
            ttk.Label(self, text=f"R{r+1}", width=4).grid(row=r+1, column=0, padx=1, pady=1)
            for c in range(self.cols):
                entry = ttk.Entry(self, width=6, justify='center')
                entry.grid(row=r+1, column=c+1, padx=1, pady=1)
                self.cells[(r, c)] = entry
                
    def add_row(self):
        if self.rows >= 50:
            return
        self.rows += 1
        r = self.rows - 1
        ttk.Label(self, text=f"R{r+1}", width=4).grid(row=r+1, column=0, padx=1, pady=1)
        for c in range(self.cols):
            entry = ttk.Entry(self, width=6, justify='center')
            entry.grid(row=r+1, column=c+1, padx=1, pady=1)
            self.cells[(r, c)] = entry
            
    def del_row(self):
        if self.rows <= 1:
            return
        r = self.rows - 1
        for c in range(self.cols):
            self.cells[(r, c)].grid_forget()
            self.cells[(r, c)].destroy()
            del self.cells[(r, c)]
        self.rows -= 1
        
    def add_col(self):
        if self.cols >= 50:
            return
        self.cols += 1
        c = self.cols - 1
        for r in range(self.rows):
            entry = ttk.Entry(self, width=6, justify='center')
            entry.grid(row=r+1, column=c+1, padx=1, pady=1)
            self.cells[(r, c)] = entry
            
    def del_col(self):
        if self.cols <= 1:
            return
        c = self.cols - 1
        for r in range(self.rows):
            self.cells[(r, c)].grid_forget()
            self.cells[(r, c)].destroy()
            del self.cells[(r, c)]
        self.cols -= 1
                
    def get_data(self):
        data = []
        for r in range(self.rows):
            row_data = []
            for c in range(self.cols):
                val = self.cells[(r, c)].get().strip()
                if val:
                    try:
                        row_data.append(float(val))
                    except ValueError:
                        row_data.append(0)
                else:
                    row_data.append(0)
            data.append(row_data)
        return np.array(data)
        
    def set_data(self, data):
        for r in range(self.rows):
            for c in range(self.cols):
                self.cells[(r, c)].delete(0, tk.END)
                
        for i, row in enumerate(data):
            for j, val in enumerate(row):
                if i < self.rows and j < self.cols:
                    self.cells[(i, j)].insert(0, str(val))


def sigmoid(x, bottom, top, ec50, hill):
    return bottom + (top - bottom) / (1 + (ec50 / x) ** hill)


def inverse_sigmoid(y, bottom, top, ec50, hill):
    return ec50 * ((top - bottom) / (y - bottom) - 1) ** (1 / hill)


class ChartWindow:
    def __init__(self, title, width=8, height=6):
        self.top = tk.Toplevel()
        self.top.title(title)
        self.top.geometry(f"{width}00x{height}00")
        
        self.fig = Figure(figsize=(width, height), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.top)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.ax = self.fig.add_subplot(111)
        
        btn_frame = ttk.Frame(self.top)
        btn_frame.pack(fill=tk.X, padx=5, pady=5)
        btn = ttk.Button(btn_frame, text='Save PNG', command=self.save_png)
        btn.pack(side=tk.LEFT, padx=5)
        add_tooltip(btn, 'Save PNG')
        
        btn = ttk.Button(btn_frame, text='Save PDF', command=self.save_pdf)
        btn.pack(side=tk.LEFT, padx=5)
        add_tooltip(btn, 'Save PDF')
        
    def save_png(self):
        filename = filedialog.asksaveasfilename(defaultextension='.png',
                                               filetypes=[('PNG图片', '*.png')])
        if filename:
            self.fig.savefig(filename, dpi=300, bbox_inches='tight')
            messagebox.showinfo('成功', f'已保存到 {filename}')
            
    def save_pdf(self):
        filename = filedialog.asksaveasfilename(defaultextension='.pdf',
                                               filetypes=[('PDF文档', '*.pdf')])
        if filename:
            self.fig.savefig(filename, dpi=300, bbox_inches='tight')
            messagebox.showinfo('成功', f'已保存到 {filename}')
        
    def clear(self):
        self.ax.clear()
        
    def show(self):
        self.fig.tight_layout()
        self.canvas.draw()


class PharmaStatApp:
    def __init__(self, root):
        self.root = root
        self.root.title('PharmaStat - Statistical Visualization')
        self.root.geometry('1200x800')
        
        self.current_indicator = tk.StringVar(value='血压(mmHg)')
        self.dose_indicator = tk.StringVar(value='抑制率(%)')
        self.corr_xlabel = tk.StringVar(value='X变量')
        self.corr_ylabel = tk.StringVar(value='Y变量')
        self.row_labels_var = tk.StringVar(value='')
        self.col_labels_var = tk.StringVar(value='')
        self.volcano_fc_var = tk.StringVar(value='2')
        self.volcano_pval_var = tk.StringVar(value='0.05')
        
        self.setup_ui()
        
        self.root.after(500, self.load_sample_data)
        
    def load_sample_data(self):
        samples = {
            'Control': [120, 122, 118, 125, 121, 119, 123, 124],
            'Model': [145, 148, 142, 150, 146, 144, 147, 149],
            'Low Dose': [138, 140, 135, 142, 139, 136, 141, 137],
            'Mid Dose': [130, 128, 132, 125, 129, 131, 127, 126],
            'High Dose': [115, 118, 112, 120, 116, 114, 117, 119],
            'Positive': [125, 128, 122, 126, 124, 127, 123, 121]
        }
        self.data_grid.set_data(samples)
        
        dose_doses = [0.001, 0.01, 0.1, 1, 10, 100, 1000]
        dose_effects = [5, 12, 35, 68, 88, 97, 100]
        dose_sds = [2, 3, 4, 5, 3, 2, 1]
        dose_ns = [6, 6, 6, 6, 6, 6, 6]
        for i, (d, e, s, n) in enumerate(zip(dose_doses, dose_effects, dose_sds, dose_ns)):
            self.dose_grid.cells[(i, 0)].insert(0, str(d))
            self.dose_grid.cells[(i, 1)].insert(0, str(e))
            self.dose_grid.cells[(i, 2)].insert(0, str(s))
            self.dose_grid.cells[(i, 3)].insert(0, str(n))
        
        np.random.seed(42)
        corr_x = np.linspace(0, 100, 30)
        corr_y = 0.5 * corr_x + 10 + np.random.normal(0, 10, 30)
        for i in range(30):
            self.corr_grid.cells[(i, 0)].insert(0, f'{corr_x[i]:.1f}')
            self.corr_grid.cells[(i, 1)].insert(0, f'{corr_y[i]:.1f}')
        
        heatmap_data = np.random.randn(10, 10)
        heatmap_data = np.cumsum(heatmap_data, axis=1)
        for r in range(10):
            for c in range(10):
                self.heatmap_grid.cells[(r, c)].insert(0, f'{heatmap_data[r,c]:.1f}')
        
        volcano_genes = [
            ('GAPDH', 0.12, 0.85),
            ('ACTB', 0.05, 0.72),
            ('TNF', 2.5, 0.001),
            ('IL6', 3.1, 0.0003),
            ('IL1B', 2.8, 0.0008),
            ('CXCL8', 1.9, 0.015),
            ('CCL2', -2.3, 0.002),
            ('CXCL10', -3.5, 0.0001),
            ('IFNG', -1.8, 0.008),
            ('IL10', -2.1, 0.004),
            ('VEGFA', 1.5, 0.025),
            ('EGFR', 0.8, 0.15),
            ('TP53', -0.3, 0.65),
            ('BAX', 0.2, 0.78),
            ('BCL2', -0.15, 0.88),
            ('CASP3', 0.4, 0.45),
            ('MMP9', 2.2, 0.003),
            ('MMP2', 1.6, 0.02),
            ('COX2', 2.9, 0.0005),
            ('NOS2', -1.5, 0.03),
        ]
        for i, (gene, fc, pval) in enumerate(volcano_genes):
            self.volcano_grid.cells[(i, 0)].insert(0, gene)
            self.volcano_grid.cells[(i, 1)].insert(0, str(fc))
            self.volcano_grid.cells[(i, 2)].insert(0, str(pval))
        
    def setup_ui(self):
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.tab_group = ttk.Frame(self.notebook)
        self.tab_dose = ttk.Frame(self.notebook)
        self.tab_corr = ttk.Frame(self.notebook)
        self.tab_heatmap = ttk.Frame(self.notebook)
        self.tab_volcano = ttk.Frame(self.notebook)
        
        self.notebook.add(self.tab_group, text='Group')
        self.notebook.add(self.tab_dose, text='Dose-Response')
        self.notebook.add(self.tab_corr, text='Correlation')
        self.notebook.add(self.tab_heatmap, text='Heatmap')
        self.notebook.add(self.tab_volcano, text='Volcano')
        
        self.setup_group_tab()
        self.setup_dose_tab()
        self.setup_corr_tab()
        self.setup_heatmap_tab()
        self.setup_volcano_tab()
        
    def setup_group_tab(self):
        input_frame = ttk.Frame(self.tab_group)
        input_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        top_bar = ttk.Frame(input_frame)
        top_bar.pack(fill=tk.X, pady=(0,5))
        
        ttk.Label(top_bar, text='Indicator:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        indicator_combo = ttk.Combobox(top_bar, textvariable=self.current_indicator, 
                                        values=DEFAULT_INDICATORS, width=18)
        indicator_combo.pack(side=tk.LEFT, padx=5)
        
        btn = ttk.Button(top_bar, text='Clear', command=self.clear_group_data)
        btn.pack(side=tk.LEFT, padx=10)
        add_tooltip(btn, 'Clear')
        
        canvas_frame = ttk.Frame(input_frame, relief=tk.SUNKEN, borderwidth=1)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        self.canvas_group = tk.Canvas(canvas_frame)
        self.vsb_group = ttk.Scrollbar(canvas_frame, orient=tk.VERTICAL, command=self.canvas_group.yview)
        self.hsb_group = ttk.Scrollbar(canvas_frame, orient=tk.HORIZONTAL, command=self.canvas_group.xview)
        self.canvas_group.configure(yscrollcommand=self.vsb_group.set, xscrollcommand=self.hsb_group.set)
        
        self.vsb_group.pack(side=tk.RIGHT, fill=tk.Y)
        self.hsb_group.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas_group.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.data_grid_frame = ttk.Frame(self.canvas_group)
        self.canvas_group.create_window((0, 0), window=self.data_grid_frame, anchor=tk.NW)
        
        self.data_grid = DataGrid(self.data_grid_frame, rows=6, cols=8)
        self.data_grid.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.data_grid_frame.bind("<Configure>", lambda e: self.canvas_group.configure(scrollregion=self.canvas_group.bbox("all")))
        
        btn_bar = ttk.Frame(input_frame)
        btn_bar.pack(fill=tk.X, pady=(5,0))
        
        btn = ttk.Button(btn_bar, text='+ Add Col', command=self.data_grid.add_col)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '+ Add Col')
        
        btn = ttk.Button(btn_bar, text='- Del Col', command=self.data_grid.del_col)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '- Del Col')
        
        ttk.Separator(btn_bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        btn = ttk.Button(btn_bar, text='Bar+Error', command=self.generate_bar_chart)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, 'Bar+Error')
        
        btn = ttk.Button(btn_bar, text='Boxplot', command=self.generate_box_chart)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, 'Boxplot')
        
        btn = ttk.Button(btn_bar, text='Import CSV', command=self.load_csv)
        btn.pack(side=tk.LEFT, padx=10)
        add_tooltip(btn, 'Import CSV')
        
        btn = ttk.Button(btn_bar, text='Export', command=self.export_stats)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, 'Export')
        
    def setup_dose_tab(self):
        left_frame = ttk.Frame(self.tab_dose)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        top_bar = ttk.Frame(left_frame)
        top_bar.pack(fill=tk.X, pady=(0,5))
        
        ttk.Label(top_bar, text='Indicator:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        ttk.Entry(top_bar, textvariable=self.dose_indicator, width=15).pack(side=tk.LEFT, padx=5)
        
        btn = ttk.Button(top_bar, text='Clear', command=self.clear_dose_data)
        btn.pack(side=tk.LEFT, padx=10)
        add_tooltip(btn, 'Clear')
        
        canvas_frame = ttk.Frame(left_frame, relief=tk.SUNKEN, borderwidth=1)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        self.canvas_dose = tk.Canvas(canvas_frame)
        self.vsb_dose = ttk.Scrollbar(canvas_frame, orient=tk.VERTICAL, command=self.canvas_dose.yview)
        self.canvas_dose.configure(yscrollcommand=self.vsb_dose.set)
        
        self.vsb_dose.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas_dose.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.dose_grid_frame = ttk.Frame(self.canvas_dose)
        self.canvas_dose.create_window((0, 0), window=self.dose_grid_frame, anchor=tk.NW)
        
        self.dose_grid = DoseResponseGrid(self.dose_grid_frame, rows=8)
        self.dose_grid.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.dose_grid_frame.bind("<Configure>", lambda e: self.canvas_dose.configure(scrollregion=self.canvas_dose.bbox("all")))
        
        btn_bar = ttk.Frame(left_frame)
        btn_bar.pack(fill=tk.X, pady=(5,0))
        
        btn = ttk.Button(btn_bar, text='+ Add Row', command=self.dose_grid.add_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '+ Add Row')
        
        btn = ttk.Button(btn_bar, text='- Del Row', command=self.dose_grid.del_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '- Del Row')
        ttk.Separator(btn_bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        btn = ttk.Button(btn_bar, text='Generate Curve', command=self.generate_dose_curve)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, 'Generate Curve')
        
    def setup_corr_tab(self):
        left_frame = ttk.Frame(self.tab_corr)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        top_bar = ttk.Frame(left_frame)
        top_bar.pack(fill=tk.X, pady=(0,5))
        
        ttk.Label(top_bar, text='X:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        ttk.Entry(top_bar, textvariable=self.corr_xlabel, width=10).pack(side=tk.LEFT, padx=5)
        ttk.Label(top_bar, text='Y:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        ttk.Entry(top_bar, textvariable=self.corr_ylabel, width=10).pack(side=tk.LEFT, padx=5)
        btn = ttk.Button(top_bar, text='Clear', command=self.clear_corr_data)
        btn.pack(side=tk.LEFT, padx=10)
        add_tooltip(btn, 'Clear')
        
        canvas_frame = ttk.Frame(left_frame, relief=tk.SUNKEN, borderwidth=1)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        self.canvas_corr = tk.Canvas(canvas_frame)
        self.vsb_corr = ttk.Scrollbar(canvas_frame, orient=tk.VERTICAL, command=self.canvas_corr.yview)
        self.canvas_corr.configure(yscrollcommand=self.vsb_corr.set)
        
        self.vsb_corr.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas_corr.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.corr_grid_frame = ttk.Frame(self.canvas_corr)
        self.canvas_corr.create_window((0, 0), window=self.corr_grid_frame, anchor=tk.NW)
        
        self.corr_grid = CorrelationGrid(self.corr_grid_frame, max_rows=30)
        self.corr_grid.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.corr_grid_frame.bind("<Configure>", lambda e: self.canvas_corr.configure(scrollregion=self.canvas_corr.bbox("all")))
        
        btn_bar = ttk.Frame(left_frame)
        btn_bar.pack(fill=tk.X, pady=(5,0))
        
        btn = ttk.Button(btn_bar, text='+ Add Row', command=self.corr_grid.add_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '+ Add Row')
        
        btn = ttk.Button(btn_bar, text='- Del Row', command=self.corr_grid.del_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '- Del Row')
        ttk.Separator(btn_bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        btn = ttk.Button(btn_bar, text='Scatter Plot', command=self.generate_corr_plot)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, 'Scatter Plot')
        
    def setup_heatmap_tab(self):
        left_frame = ttk.Frame(self.tab_heatmap)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        top_bar = ttk.Frame(left_frame)
        top_bar.pack(fill=tk.X, pady=(0,5))
        
        ttk.Label(top_bar, text='Row Labels:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        ttk.Entry(top_bar, textvariable=self.row_labels_var, width=20).pack(side=tk.LEFT, padx=5)
        ttk.Label(top_bar, text='Col Labels:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        ttk.Entry(top_bar, textvariable=self.col_labels_var, width=20).pack(side=tk.LEFT, padx=5)
        btn = ttk.Button(top_bar, text='Clear', command=self.clear_heatmap_data)
        btn.pack(side=tk.LEFT, padx=10)
        add_tooltip(btn, 'Clear')
        
        canvas_frame = ttk.Frame(left_frame, relief=tk.SUNKEN, borderwidth=1)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        self.canvas_heatmap = tk.Canvas(canvas_frame)
        self.vsb_heatmap = ttk.Scrollbar(canvas_frame, orient=tk.VERTICAL, command=self.canvas_heatmap.yview)
        self.hsb_heatmap = ttk.Scrollbar(canvas_frame, orient=tk.HORIZONTAL, command=self.canvas_heatmap.xview)
        self.canvas_heatmap.configure(yscrollcommand=self.vsb_heatmap.set, xscrollcommand=self.hsb_heatmap.set)
        
        self.vsb_heatmap.pack(side=tk.RIGHT, fill=tk.Y)
        self.hsb_heatmap.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas_heatmap.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.heatmap_grid_frame = ttk.Frame(self.canvas_heatmap)
        self.canvas_heatmap.create_window((0, 0), window=self.heatmap_grid_frame, anchor=tk.NW)
        
        self.heatmap_grid = HeatmapGrid(self.heatmap_grid_frame, rows=10, cols=10)
        self.heatmap_grid.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.heatmap_grid_frame.bind("<Configure>", lambda e: self.canvas_heatmap.configure(scrollregion=self.canvas_heatmap.bbox("all")))
        
        btn_bar = ttk.Frame(left_frame)
        btn_bar.pack(fill=tk.X, pady=(5,0))
        
        btn = ttk.Button(btn_bar, text='+ Add Row', command=self.heatmap_grid.add_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '+ Add Row')
        
        btn = ttk.Button(btn_bar, text='- Del Row', command=self.heatmap_grid.del_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '- Del Row')
        
        btn = ttk.Button(btn_bar, text='+ Add Col', command=self.heatmap_grid.add_col)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '+ Add Col')
        
        btn = ttk.Button(btn_bar, text='- Del Col', command=self.heatmap_grid.del_col)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '- Del Col')
        ttk.Separator(btn_bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        btn = ttk.Button(btn_bar, text='Generate', command=self.generate_heatmap)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, 'Generate')
        
    def setup_volcano_tab(self):
        left_frame = ttk.Frame(self.tab_volcano)
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        top_bar = ttk.Frame(left_frame)
        top_bar.pack(fill=tk.X, pady=(0,5))
        
        ttk.Label(top_bar, text='FC Threshold:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        ttk.Entry(top_bar, textvariable=self.volcano_fc_var, width=10).pack(side=tk.LEFT, padx=5)
        
        ttk.Label(top_bar, text='P-value:', font=('Arial', 11)).pack(side=tk.LEFT, padx=5)
        ttk.Entry(top_bar, textvariable=self.volcano_pval_var, width=10).pack(side=tk.LEFT, padx=5)
        
        btn = ttk.Button(top_bar, text='Clear', command=self.clear_volcano_data)
        btn.pack(side=tk.LEFT, padx=10)
        add_tooltip(btn, 'Clear')
        
        canvas_frame = ttk.Frame(left_frame, relief=tk.SUNKEN, borderwidth=1)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        
        self.canvas_volcano = tk.Canvas(canvas_frame)
        self.vsb_volcano = ttk.Scrollbar(canvas_frame, orient=tk.VERTICAL, command=self.canvas_volcano.yview)
        self.canvas_volcano.configure(yscrollcommand=self.vsb_volcano.set)
        
        self.vsb_volcano.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas_volcano.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        self.volcano_grid_frame = ttk.Frame(self.canvas_volcano)
        self.canvas_volcano.create_window((0, 0), window=self.volcano_grid_frame, anchor=tk.NW)
        
        self.volcano_grid = VolcanoGrid(self.volcano_grid_frame, rows=20)
        self.volcano_grid.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        self.volcano_grid_frame.bind("<Configure>", lambda e: self.canvas_volcano.configure(scrollregion=self.canvas_volcano.bbox("all")))
        
        btn_bar = ttk.Frame(left_frame)
        btn_bar.pack(fill=tk.X, pady=(5,0))
        
        btn = ttk.Button(btn_bar, text='+ Add Row', command=self.volcano_grid.add_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '+ Add Row')
        
        btn = ttk.Button(btn_bar, text='- Del Row', command=self.volcano_grid.del_row)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, '- Del Row')
        
        ttk.Separator(btn_bar, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        
        btn = ttk.Button(btn_bar, text='Generate Volcano', command=self.generate_volcano_plot)
        btn.pack(side=tk.LEFT, padx=3)
        add_tooltip(btn, 'Generate Volcano')
        
    def clear_volcano_data(self):
        for r in range(self.volcano_grid.rows):
            for c in range(3):
                self.volcano_grid.cells[(r, c)].delete(0, tk.END)
        
    def clear_group_data(self):
        for r in range(6):
            for c in range(self.data_grid.cols):
                self.data_grid.cells[(r, c)].delete(0, tk.END)
                
    def clear_dose_data(self):
        for r in range(self.dose_grid.rows):
            for c in range(4):
                self.dose_grid.cells[(r, c)].delete(0, tk.END)
                
    def clear_corr_data(self):
        for r in range(self.corr_grid.max_rows):
            for c in range(3):
                self.corr_grid.cells[(r, c)].delete(0, tk.END)
                
    def clear_heatmap_data(self):
        for r in range(self.heatmap_grid.rows):
            for c in range(self.heatmap_grid.cols):
                self.heatmap_grid.cells[(r, c)].delete(0, tk.END)
                
    def collect_data(self):
        return self.data_grid.get_data()
        
    def calculate_stats(self, data):
        stats_dict = {}
        for group, values in data.items():
            values = np.array(values)
            stats_dict[group] = {
                'mean': np.mean(values),
                'std': np.std(values, ddof=1),
                'sem': np.std(values, ddof=1) / np.sqrt(len(values)),
                'n': len(values),
                'median': np.median(values),
            }
        return stats_dict
        
    def statistical_test(self, data, group):
        if group == 'Control' or 'Control' not in data:
            return None
        if len(data[group]) < 2:
            return None
        t_stat, p_value = stats.ttest_ind(data[group], data['Control'])
        return p_value
        
    def get_significance_marker(self, p_value):
        if p_value is None:
            return ''
        if p_value < 0.001:
            return '***'
        elif p_value < 0.01:
            return '**'
        elif p_value < 0.05:
            return '*'
        return ''
        
    def generate_bar_chart(self):
        data = self.collect_data()
        if not data:
            messagebox.showwarning('警告', '请输入数据')
            return
            
        stats_dict = self.calculate_stats(data)
        
        win = ChartWindow(f'{self.current_indicator.get()} - 柱状图')
        ax = win.ax
        
        groups = [g for g in GROUP_ORDER if g in data]
        x_pos = np.arange(len(groups))
        means = [stats_dict[g]['mean'] for g in groups]
        sems = [stats_dict[g]['sem'] for g in groups]
        colors = [GROUP_COLORS[GROUP_ORDER.index(g)] for g in groups]
        
        bars = ax.bar(x_pos, means, yerr=sems, capsize=5, color=colors, 
                          edgecolor='black', linewidth=1.2, error_kw={'linewidth': 1.5})
        
        for i, group in enumerate(groups):
            p_value = self.statistical_test(data, group)
            marker = self.get_significance_marker(p_value)
            if marker and group != 'Control':
                max_val = stats_dict[group]['mean'] + stats_dict[group]['sem']
                ax.annotate(marker, (x_pos[i], max_val + sems[i]*0.3),
                               ha='center', va='bottom', fontsize=14, fontweight='bold')
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(groups, fontsize=11, rotation=15)
        ax.set_ylabel(self.current_indicator.get(), fontsize=12)
        ax.set_title(f'{self.current_indicator.get()} - 柱状图 (均值±SEM)', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        
        win.show()
        
    def generate_box_chart(self):
        data = self.collect_data()
        if not data:
            messagebox.showwarning('警告', '请输入数据')
            return
            
        win = ChartWindow(f'{self.current_indicator.get()} - 箱线图')
        ax = win.ax
        
        groups = [g for g in GROUP_ORDER if g in data]
        plot_data = [data[g] for g in groups]
        colors = [GROUP_COLORS[GROUP_ORDER.index(g)] for g in groups]
        
        bp = ax.boxplot(plot_data, positions=range(len(groups)), 
                           widths=0.6, patch_artist=True)
        
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
            
        for i, group in enumerate(groups):
            p_value = self.statistical_test(data, group)
            marker = self.get_significance_marker(p_value)
            if marker and group != 'Control':
                y_max = max(data[group])
                ax.annotate(marker, (i, y_max + (y_max*0.05)),
                               ha='center', va='bottom', fontsize=14, fontweight='bold')
        
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, fontsize=11, rotation=15)
        ax.set_ylabel(self.current_indicator.get(), fontsize=12)
        ax.set_title(f'{self.current_indicator.get()} - 箱线图', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        
        win.show()
        
    def generate_dose_curve(self):
        doses, effects, sds, ns = self.dose_grid.get_data()
        if len(doses) < 4:
            messagebox.showwarning('警告', '请输入至少4组数据')
            return
            
        win = ChartWindow('量效曲线')
        ax = win.ax
        
        doses = np.array(doses)
        effects = np.array(effects)
        
        x_fit = np.logspace(np.log10(min(doses)), np.log10(max(doses)), 100)
        
        try:
            p0 = [min(effects), max(effects), np.median(doses), 2]
            popt, _ = curve_fit(sigmoid, doses, effects, p0=p0, maxfev=5000)
            
            bottom, top, ec50, hill = popt
            y_fit = sigmoid(x_fit, *popt)
            
            ax.plot(x_fit, y_fit, 'b-', linewidth=2, label='Sigmoid fit')
            
            if bottom < 50 < top:
                ic50 = inverse_sigmoid(50, *popt)
                ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5)
                ax.axvline(x=ic50, color='red', linestyle='--', alpha=0.5)
                ax.annotate(f'IC50={ic50:.2f}', (ic50, 50), xytext=(5, 5), textcoords='offset points', fontsize=10)
            
            info_text = f'IC50={ec50:.2f}\nHill={hill:.2f}'
            ax.text(0.05, 0.95, info_text, transform=ax.transAxes, fontsize=10,
                       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        except:
            ax.plot(doses, effects, 'b-', linewidth=2, label='Data')
        
        ax.errorbar(doses, effects, yerr=sds, fmt='o', capsize=5, color='red', 
                          ecolor='gray', markersize=8, label='Data')
        
        ax.set_xscale('log')
        ax.set_xlabel('Concentration / Dose', fontsize=12)
        ax.set_ylabel(self.dose_indicator.get(), fontsize=12)
        ax.set_title('Dose-Response Curve', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.legend()
        
        win.show()
        
    def generate_corr_plot(self):
        x_data, y_data, groups = self.corr_grid.get_data()
        if len(x_data) < 3:
            messagebox.showwarning('警告', '请输入至少3组数据')
            return
            
        win = ChartWindow('相关性散点图')
        ax = win.ax
        
        x_data = np.array(x_data)
        y_data = np.array(y_data)
        
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_data, y_data)
        
        x_line = np.linspace(min(x_data), max(x_data), 100)
        y_line = slope * x_line + intercept
        
        ax.scatter(x_data, y_data, c='steelblue', s=60, alpha=0.7, edgecolors='white')
        ax.plot(x_line, y_line, 'r-', linewidth=2, label=f'y={slope:.2f}x+{intercept:.2f}')
        
        text = f'r={r_value:.3f}\np={p_value:.4f}'
        if p_value < 0.001:
            text += '***'
        elif p_value < 0.01:
            text += '**'
        elif p_value < 0.05:
            text += '*'
            
        ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=11,
                         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_xlabel(self.corr_xlabel.get(), fontsize=12)
        ax.set_ylabel(self.corr_ylabel.get(), fontsize=12)
        ax.set_title('Correlation Plot', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.legend()
        
        win.show()
        
    def generate_heatmap(self):
        data = self.heatmap_grid.get_data()
        
        if np.all(data == 0):
            messagebox.showwarning('警告', '请输入数据')
            return
            
        win = ChartWindow('热图')
        ax = win.ax
        
        row_labels = self.row_labels_var.get().split(',') if self.row_labels_var.get() else [f'R{i+1}' for i in range(data.shape[0])]
        col_labels = self.col_labels_var.get().split(',') if self.col_labels_var.get() else [f'C{i+1}' for i in range(data.shape[1])]
        
        if row_labels and len(row_labels) != data.shape[0]:
            row_labels = [f'R{i+1}' for i in range(data.shape[0])]
        if col_labels and len(col_labels) != data.shape[1]:
            col_labels = [f'C{i+1}' for i in range(data.shape[1])]
        
        df_heatmap = pd.DataFrame(data, index=row_labels, columns=col_labels)
        
        sns.heatmap(df_heatmap, annot=True, fmt='.1f', cmap='coolwarm', 
                    center=0, linewidths=1, linecolor='white',
                    ax=ax, cbar_kws={'shrink': 0.8})
        
        ax.set_title('Heatmap', fontsize=14, fontweight='bold')
        
        win.show()
        
    def generate_volcano_plot(self):
        genes, log2fcs, pvals = self.volcano_grid.get_data()
        if len(genes) < 3:
            messagebox.showwarning('Warning', 'Please input at least 3 genes')
            return
            
        log2fcs = np.array(log2fcs)
        pvals = np.array(pvals)
        neg_log10_pvals = -np.log10(pvals + 1e-10)
        
        fc_threshold = self.volcano_fc_var.get()
        pval_threshold = self.volcano_pval_var.get()
        
        colors = ['gray'] * len(genes)
        for i in range(len(genes)):
            if log2fcs[i] > np.log2(float(fc_threshold)) and pvals[i] < float(pval_threshold):
                colors[i] = '#E53935'
            elif log2fcs[i] < -np.log2(float(fc_threshold)) and pvals[i] < float(pval_threshold):
                colors[i] = '#1976D2'
        
        win = ChartWindow('Volcano Plot')
        ax = win.ax
        
        ax.scatter(log2fcs, neg_log10_pvals, c=colors, alpha=0.7, s=50, edgecolors='white')
        
        ax.axhline(y=-np.log10(float(pval_threshold)), color='blue', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(x=np.log2(float(fc_threshold)), color='red', linestyle='--', alpha=0.5, linewidth=1)
        ax.axvline(x=-np.log2(float(fc_threshold)), color='red', linestyle='--', alpha=0.5, linewidth=1)
        
        for i, gene in enumerate(genes):
            if abs(log2fcs[i]) > np.log2(float(fc_threshold)) and pvals[i] < float(pval_threshold):
                ax.annotate(gene, (log2fcs[i], neg_log10_pvals[i]), fontsize=8, alpha=0.8)
        
        ax.set_xlabel('Log2(Fold Change)', fontsize=12)
        ax.set_ylabel('-Log10(P-value)', fontsize=12)
        ax.set_title('Volcano Plot', fontsize=14, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        win.show()
        
    def load_csv(self):
        filename = filedialog.askopenfilename(filetypes=[('CSV文件', '*.csv'), ('所有文件', '*.*')])
        if filename:
            try:
                df = pd.read_csv(filename)
                data = {}
                for group in GROUP_ORDER:
                    if group in df.columns:
                        data[group] = df[group].dropna().tolist()
                    else:
                        data[group] = []
                self.data_grid.set_data(data)
                if len(df.columns) > 0:
                    self.current_indicator.set(df.columns[0])
                messagebox.showinfo('成功', f'已导入 {filename}')
            except Exception as e:
                messagebox.showerror('错误', str(e))
                
    def export_stats(self):
        data = self.collect_data()
        if not data:
            messagebox.showwarning('警告', '请输入数据')
            return
            
        stats_dict = self.calculate_stats(data)
        
        rows = []
        for group in GROUP_ORDER:
            if group in stats_dict:
                s = stats_dict[group]
                rows.append({
                    'Group': group,
                    'N': s['n'],
                    'Mean': f"{s['mean']:.2f}",
                    'SD': f"{s['std']:.2f}",
                    'SEM': f"{s['sem']:.2f}",
                })
        
        df = pd.DataFrame(rows)
        
        filename = filedialog.asksaveasfilename(defaultextension='.csv',
                                               filetypes=[('CSV文件', '*.csv')])
        if filename:
            df.to_csv(filename, index=False, encoding='utf-8-sig')
            messagebox.showinfo('成功', f'已导出到 {filename}')


def main():
    root = tk.Tk()
    app = PharmaStatApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()