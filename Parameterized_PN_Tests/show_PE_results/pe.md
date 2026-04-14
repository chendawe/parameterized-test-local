Truth: h(t; $\theta^\text{target}$)，实际上对应信号数据

Template: h(t; $\theta$)
1. 从什么分布里面取？
2. 怎么猜？
3. 怎么衡量相似度？


$\{m_1, m_2; \chi_{1L}， \chi_{2L}； \dots\}$
先验分布（prior）


$f(t_n), f(t_{n+1})$

MCMC(Markov Chain Monte Carlo)
采样，MCMC（采样空间：初始点，下一步->怎么走？->每一步往哪个方向走多少->趋向目标波形参数(高似然区域)->完成收敛(基于终止条件)）

似然：(1,0),(0.5,1)，归一化之后做内积；
噪声noise fd PSD，权重 $\propto 1/\text{PSD}(f)$
h(t) -> $\tilde{h}$(f)

PSD(Power Spectral Density)

$s(t) = h(t) + n(t)$

$s(f) = h(f) + n(f)$

$\braket{h(t)|g(t)}=\braket{\tilde{h}(f)|\tilde{g}(f)}
= \int \left[ \frac{\tilde{h}(f)*\tilde{g}(f)}{\text{PSD}(f)} \right]\text{d}f$

带入高斯噪声模型：
似然：Likelihood $\mathcal{L} \propto e^{\text{内积}}$


=================


PE parameter estimation

$\{m_1, m_2; \chi_{1L}， \chi_{2L}； \dots\}$里面每个参数的分布结果（后验分布posterior），带有不确定性的，要进行统计学分析的分布


==============


线性回归：
f(x) = kx + b



1. 现象学频域FD波形的构造(蕴含GR信息) -> 依赖$\theta$的GR波形模板
2. GR模板 -> 参数估计 -> $\theta$的后验分布(GR预测结果)
3. 参数化检验GR

$h(f; \theta) = h(f; \theta^{(\text{GR})})$

IMRPhenomXAS，AS = Aligned Spin
- $h(f; \theta^{(\text{GR})}) = h(f; \theta^{(\text{GR})}, \psi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L}))$
    - $\theta^{(\text{GR})} = \{m_1, m_2; \chi_{1L}， \chi_{2L}； \dots\}$
- 参数化：$h(f; \theta^{(\text{Parameterized})}) = h(f; \psi_{j(l)}\in\theta^{(\text{Parameterized})})$
    - $\theta^{(\text{Parameterized})} = \{m_1, m_2; \chi_{1L}， \chi_{2L}；\psi_{j(l)}; \dots\}$


一堆二维数据点
- 有个人的理论说： $a = a(b,c) = b^2-c^2$
    - $f(x) = (b^2-c^2)\cdot x^2 + b\cdot x + c$
- 那这样检验：$f(x) = a\cdot x^2 + b\cdot x + c$
    - 检查：$a = a(b,c) = b^2-c^2$?

已有的引力波信号(strain)，gwf数据文件
- GR: $\psi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L})$
    - $h(f; \theta^{(\text{GR})})$
- $h(f;\psi_{j(l)}\in\theta^{(\text{Parameterized})})$
    - $\psi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L})$，GR公式算出来的
    - $\psi_{j(l)}$，作为独立参数估计出来的



================

LIGO：
- fractional: $\phi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L}) \longrightarrow \phi(m_1,m_2,\chi_{1L},\chi_{2L}) \cdot (1 + \delta_{j(l)})$
    - $\phi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L}) + \Delta_{j(l)} = \phi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L}) \cdot (1 + \frac{\Delta_{j(l)}}{\phi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L})})$

我们：
- Absoslute：$\psi_{j(l)}(m_1,m_2,\chi_{1L},\chi_{2L}) \longrightarrow \psi_{j(l)}$


============

0,2,6l
- $\psi_{j(l)}(m_1,m_2) = \psi_{j(l)}$
    - if 固定 $m_1$，变动$m_2$，根据GR的约束，$\psi_{j(l)}$随$m_2$以GR规则变动；

- 反解$m_2 \mapsto \psi_{j(l)}$，得到$\psi_{j(l)} \mapsto m_2$

- 在$m_1$-$m_2$平面内
    - 已有$m_1$,$m_2$的样本，可以画出可信区域


=========

可信区域
- 对称式
- 概率密度
    - 二维样本点，可以在平面上画出联合分布（散点）
    - 核密度估计KDE来估计得到连续的联合概率密度分布
    - 概率密度的等高线画出闭合边界，闭合边界围成可信区域，占总样本数后验概率$p$的样本点落在可信区域内，这个可信区域对应后验概率$p$


======

- 估计到的$m_2$算的$\psi$，和估计的$\psi$一样；
- 估计到的$\psi$算的$m_2$，和估计的$m_2$一样；

========

- 直接画两个质量的样本点：$\{(m_1,m_2)\}$
- 表达$\psi$：$\{(m_1,m_2(\psi))\}$

0.4，{0.3,0.5,0.7},0.7-0.3 = 0.4

s1 = {m1,m2,chi1,chi2}
s2 = {m1,m2,chi1,chi2}
s3 = {m1,m2,chi1,chi2}
...




- 用不同替代理论的关系式
- 量化框架优化：
    - 对自旋依赖的处理
- 不一定是PN系数，波形中其它参数结构的参数化
- 可信区间的构造，先验分布的选择
- SNR影响


4-OGC，IMRPhenomXPHM，ini文件