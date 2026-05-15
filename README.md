<div align="center">

# 🐉 板凳龙闹元宵

### 2024年高教社杯全国大学生数学建模竞赛 A题 完整求解

[![Python](https://img.shields.io/badge/Python-3.8+-blue?logo=python&logoColor=white)](https://www.python.org)
[![NumPy](https://img.shields.io/badge/NumPy-013243?logo=numpy&logoColor=white)](https://numpy.org)
[![SciPy](https://img.shields.io/badge/SciPy-8CAAE6?logo=scipy&logoColor=white)](https://scipy.org)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![Version](https://img.shields.io/badge/Version-2.4-orange)](https://github.com)
[![CUMCM](https://img.shields.io/badge/CUMCM-2024A-red)](https://github.com)

**🤖 由 Claude Code + DeepSeek V4 辅助生成**

</div>

---

## 📋 项目概述

本仓库包含2024年高教社杯全国大学生数学建模竞赛 **A题"板凳龙"闹元宵** 的完整求解方案。

**问题背景：** 浙闽地区传统"板凳龙"民俗活动——将多条板凳首尾相连形成龙形，龙头沿螺旋线盘入盘出。需要模拟龙队运动、检测碰撞、优化螺距、设计调头曲线并确定最大安全速度。

### 求解的五个问题

| 问题 | 描述 | 方法 |
|:----:|------|------|
| **1** | 基础盘入模拟（0~300s，螺距55cm） | 链式约束求解 + 数值弧长反算 |
| **2** | 碰撞检测与终止时刻 | 多边形相交检测 + AABB预检 |
| **3** | 最小螺距确定 | 二分法搜索 + 可行性验证 |
| **4** | S形调头曲线设计（-100s~100s） | 双圆弧调头路径 + 盘出螺线中心对称 |
| **5** | 最大龙头速度确定 | 速度放大因子分析 |

---

## 🏗️ 代码架构

```
├── bench_dragon_simulation.py   # 主求解代码 (v2.4优化版)
├── bench_dragon_comprehensive.md # 完整综合文档
├── problem_analysis.md           # 题目分析文档
└── README.md                     # 本文件
```

### 代码模块

| 模块 | 说明 |
|:----|:----|
| `第0部分` | 基础参数、阿基米德螺线几何函数、链式约束求解 |
| `第1部分` | `DragonSimulation` 模拟引擎（含hot-start优化） |
| `第1.5部分` | `CompositePath` 组合路径 + 统一路径弦长约束求解器 |
| `第2部分` | 问题1：基础盘入模拟 |
| `第3部分` | 问题2：碰撞检测 + 终止时刻 |
| `第4部分` | 问题3：最小螺距二分搜索 |
| `第5部分` | 问题4：S形调头曲线 + 完整模拟 |
| `第6部分` | 问题5：最大龙头速度分析 |

---

## 🚀 快速开始

### 环境要求

```bash
# Python 3.8+
pip install numpy scipy openpyxl matplotlib
```

### 运行

```bash
python bench_dragon_simulation.py
```

### 预期输出

| 文件 | 说明 |
|:----|:----|
| `result1.xlsx` | 问题1：0~300s盘入模拟 |
| `result2.xlsx` | 问题2：终止时刻状态 |
| `result3.xlsx` | 问题3：最小螺距结果 |
| `result4.xlsx` | 问题4：-100s~100s调头模拟 |
| `result5_speed_analysis.xlsx` | 问题5：速度分析 |
| `simulation_plots.png` | 轨迹可视化 |

---

## 🔧 关键技术细节

### 阿基米德螺线

$$r = b\theta, \quad b = \frac{p}{2\pi}$$

弧长解析公式：
$$s(\theta) = \frac{b}{2}\left[\theta\sqrt{\theta^2+1} + \sinh^{-1}(\theta)\right]$$

### 链式约束（核心算法）

各把手均位于螺线上，相邻把手满足弦长约束：
$$|P(\theta_i) - P(\theta_{i+1})| = L_i$$

其中 $L_0 = 2.86\text{m}$（龙头），$L_{1..222} = 1.65\text{m}$（龙身）。

### 盘出螺线

盘出螺线是盘入螺线的**中心对称像**：
$$E(\theta) = -P(\theta) \quad \text{(注意：} \neq P(\theta+\pi)\text{，相差 }b\pi \approx 0.85\text{m/圈)}$$

---

## 📜 版本历史

| 版本 | 要点 |
|:----|:----|
| **v2.4** | 八大性能优化：hot-start链求解、theta直接求解、紧区间、全覆盖量化等 |
| **v2.3** | 盘出螺线中心对称修正、S曲线C¹连续性再修正 |
| **v2.2** | 碰撞检测CW/CCW自适应、边交叉检测 |
| **v2.1** | 弦长约束修正、碰撞检测全范围扩展 |
| **v2.0** | 初始完整求解框架 |

---

## 🤖 AI 辅助说明

本项目代码由 **Claude Code** 接入 **DeepSeek V4** 模型辅助生成，经过多轮审查、优化和修正，最终形成完整的数学建模求解方案。

---

<div align="center">

**⭐ 如果这个项目有帮助，欢迎Star！**

📅 2026年5月

</div>
