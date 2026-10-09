<div class="sc-hero">
<p class="sc-kicker">Gusev Lab · Dana-Farber Cancer Institute</p>
<h1 class="sc-title">sc<span>HINT</span></h1>
<p class="sc-tagline">Single-cell gene expression <strong>H</strong>eritability with gene-by-cell-context <strong>INT</strong>eractions</p>
<p class="sc-lead">Estimate how much of a gene's expression is genetic, and how much of that genetic effect changes with cell state, cell type or perturbation, directly from single-cell data.</p>
<p class="sc-actions">
<a class="sc-btn sc-btn-primary" href="articles/scHINT.html">Get started</a>
<a class="sc-btn" href="articles/model.html">The model</a>
<a class="sc-btn" href="https://github.com/LidaWangPSU/scHINT">GitHub</a>
</p>
</div>


<div class="sc-figure">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 980 330" role="img" aria-label="scHINT overview: cells nested in donors, pairs of cells, Haseman-Elston regression, variance components" font-family="Nunito, 'Helvetica Neue', Arial, sans-serif">
<style>.t{font-size:13px;fill:#031634;font-weight:800;letter-spacing:.6px;text-transform:uppercase}.s{font-size:12px;fill:#51607a}.m{font-family:"Fira Mono",Menlo,monospace;font-size:11.5px;fill:#031634}.l{font-size:12px;fill:#031634;font-weight:700}</style>
<rect width="980" height="330" rx="14" fill="#f6f9fd" stroke="#d5deea"/>
<text class="t" x="34" y="38">Data</text>
<text class="t" x="300" y="38">Pairs of cells</text>
<text class="t" x="560" y="38">Regression</text>
<text class="t" x="790" y="38">Output</text>
<rect x="34" y="62" width="200" height="62" rx="10" fill="#fff" stroke="#d5deea"/>
<text class="s" x="44" y="78">donor i</text>
<rect x="44" y="105" width="4" height="9" fill="#9aa7b8"/>
<rect x="51" y="104" width="4" height="10" fill="#9aa7b8"/>
<rect x="58" y="107" width="4" height="7" fill="#9aa7b8"/>
<rect x="65" y="97" width="4" height="17" fill="#9aa7b8"/>
<rect x="72" y="102" width="4" height="12" fill="#9aa7b8"/>
<rect x="79" y="101" width="4" height="13" fill="#9aa7b8"/>
<rect x="86" y="106" width="4" height="8" fill="#9aa7b8"/>
<circle cx="110.3" cy="107.2" r="5.2" fill="#bad9f6" stroke="#fff" stroke-width="1"/>
<circle cx="114.8" cy="91.3" r="5.2" fill="#b2d4f3" stroke="#fff" stroke-width="1"/>
<circle cx="118.5" cy="110.4" r="5.2" fill="#acd0f1" stroke="#fff" stroke-width="1"/>
<circle cx="141.6" cy="84.7" r="5.2" fill="#87b6e4" stroke="#fff" stroke-width="1"/>
<circle cx="149.8" cy="101.9" r="5.2" fill="#79addf" stroke="#fff" stroke-width="1"/>
<circle cx="168.4" cy="91.3" r="5.2" fill="#5b99d4" stroke="#fff" stroke-width="1"/>
<circle cx="171.7" cy="91.6" r="5.2" fill="#5595d2" stroke="#fff" stroke-width="1"/>
<circle cx="198.5" cy="88.6" r="5.2" fill="#2978c2" stroke="#fff" stroke-width="1"/>
<circle cx="222.8" cy="92.1" r="5.2" fill="#015db4" stroke="#fff" stroke-width="1"/>
<rect x="34" y="138" width="200" height="62" rx="10" fill="#fff" stroke="#d5deea"/>
<text class="s" x="44" y="154">donor j</text>
<rect x="44" y="173" width="4" height="17" fill="#9aa7b8"/>
<rect x="51" y="179" width="4" height="11" fill="#9aa7b8"/>
<rect x="58" y="183" width="4" height="7" fill="#9aa7b8"/>
<rect x="65" y="175" width="4" height="15" fill="#9aa7b8"/>
<rect x="72" y="179" width="4" height="11" fill="#9aa7b8"/>
<rect x="79" y="174" width="4" height="16" fill="#9aa7b8"/>
<rect x="86" y="178" width="4" height="12" fill="#9aa7b8"/>
<circle cx="118.4" cy="176.0" r="5.2" fill="#add0f1" stroke="#fff" stroke-width="1"/>
<circle cx="128.6" cy="184.7" r="5.2" fill="#9cc5eb" stroke="#fff" stroke-width="1"/>
<circle cx="142.8" cy="183.7" r="5.2" fill="#84b5e3" stroke="#fff" stroke-width="1"/>
<circle cx="162.9" cy="174.2" r="5.2" fill="#649fd7" stroke="#fff" stroke-width="1"/>
<circle cx="166.7" cy="171.6" r="5.2" fill="#5d9bd5" stroke="#fff" stroke-width="1"/>
<circle cx="171.5" cy="176.8" r="5.2" fill="#5595d2" stroke="#fff" stroke-width="1"/>
<circle cx="208.4" cy="172.1" r="5.2" fill="#196dbd" stroke="#fff" stroke-width="1"/>
<circle cx="213.4" cy="164.5" r="5.2" fill="#1167ba" stroke="#fff" stroke-width="1"/>
<circle cx="216.4" cy="168.5" r="5.2" fill="#0c64b8" stroke="#fff" stroke-width="1"/>
<rect x="34" y="214" width="200" height="62" rx="10" fill="#fff" stroke="#d5deea"/>
<text class="s" x="44" y="230">donor k</text>
<rect x="44" y="248" width="4" height="18" fill="#9aa7b8"/>
<rect x="51" y="260" width="4" height="6" fill="#9aa7b8"/>
<rect x="58" y="259" width="4" height="7" fill="#9aa7b8"/>
<rect x="65" y="260" width="4" height="6" fill="#9aa7b8"/>
<rect x="72" y="253" width="4" height="13" fill="#9aa7b8"/>
<rect x="79" y="250" width="4" height="16" fill="#9aa7b8"/>
<rect x="86" y="256" width="4" height="10" fill="#9aa7b8"/>
<circle cx="115.7" cy="241.1" r="5.2" fill="#b1d3f3" stroke="#fff" stroke-width="1"/>
<circle cx="124.8" cy="248.2" r="5.2" fill="#a2c9ee" stroke="#fff" stroke-width="1"/>
<circle cx="159.2" cy="252.5" r="5.2" fill="#6aa3da" stroke="#fff" stroke-width="1"/>
<circle cx="168.2" cy="253.8" r="5.2" fill="#5b99d4" stroke="#fff" stroke-width="1"/>
<circle cx="181.7" cy="241.6" r="5.2" fill="#458acc" stroke="#fff" stroke-width="1"/>
<circle cx="183.2" cy="245.1" r="5.2" fill="#4289cb" stroke="#fff" stroke-width="1"/>
<circle cx="186.2" cy="259.5" r="5.2" fill="#3d85ca" stroke="#fff" stroke-width="1"/>
<circle cx="189.3" cy="255.8" r="5.2" fill="#3882c8" stroke="#fff" stroke-width="1"/>
<circle cx="214.0" cy="243.8" r="5.2" fill="#1067b9" stroke="#fff" stroke-width="1"/>
<text class="s" x="34" y="296">genotype g (bars) · cell state c (light to dark)</text>
<path d="M250 170 L286 170" stroke="#0074d9" stroke-width="2" marker-end="url(#ar)"/>
<defs><marker id="ar" markerWidth="9" markerHeight="9" refX="7" refY="4.5" orient="auto"><path d="M0 0 L9 4.5 L0 9 z" fill="#0074d9"/></marker></defs>
<rect x="300" y="62" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="322" y="62" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="344" y="62" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="366" y="62" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="388" y="62" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="410" y="62" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="432" y="62" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="454" y="62" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="476" y="62" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="300" y="84" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="322" y="84" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="344" y="84" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="366" y="84" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="388" y="84" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="410" y="84" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="432" y="84" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="454" y="84" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="476" y="84" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="300" y="106" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="322" y="106" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="344" y="106" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="366" y="106" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="388" y="106" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="410" y="106" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="432" y="106" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="454" y="106" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="476" y="106" width="20" height="20" rx="3" fill="#0074d9" opacity="0.20"/>
<rect x="300" y="128" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="322" y="128" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="344" y="128" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="366" y="128" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="388" y="128" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="410" y="128" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="432" y="128" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="454" y="128" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="476" y="128" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="300" y="150" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="322" y="150" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="344" y="150" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="366" y="150" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="388" y="150" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="410" y="150" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="432" y="150" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="454" y="150" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="476" y="150" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="300" y="172" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="322" y="172" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="344" y="172" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="366" y="172" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="388" y="172" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="410" y="172" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="432" y="172" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="454" y="172" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="476" y="172" width="20" height="20" rx="3" fill="#0074d9" opacity="0.37"/>
<rect x="300" y="194" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="322" y="194" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="344" y="194" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="366" y="194" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="388" y="194" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="410" y="194" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="432" y="194" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="454" y="194" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="476" y="194" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="300" y="216" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="322" y="216" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="344" y="216" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="366" y="216" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="388" y="216" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="410" y="216" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="432" y="216" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="454" y="216" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<rect x="476" y="216" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="300" y="238" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="322" y="238" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="344" y="238" width="20" height="20" rx="3" fill="#0074d9" opacity="0.46"/>
<rect x="366" y="238" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="388" y="238" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="410" y="238" width="20" height="20" rx="3" fill="#0074d9" opacity="0.29"/>
<rect x="432" y="238" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="454" y="238" width="20" height="20" rx="3" fill="#009999" opacity="0.85"/>
<rect x="476" y="238" width="20" height="20" rx="3" fill="#f6f9fd" opacity="1.00"/>
<text class="s" x="300" y="280">product of expression for every pair of cells</text>
<rect x="300" y="294" width="11" height="11" rx="2" fill="#009999" opacity=".85"/><text class="s" x="317" y="304">same donor</text>
<rect x="404" y="294" width="11" height="11" rx="2" fill="#0074d9" opacity=".4"/><text class="s" x="421" y="304">related donors</text>
<path d="M515 170 L545 170" stroke="#0074d9" stroke-width="2" marker-end="url(#ar)"/>
<text class="m" x="560" y="72">y<tspan baseline-shift="sub" font-size="9">im</tspan> y<tspan baseline-shift="sub" font-size="9">jn</tspan> ~</text>
<rect x="560" y="81" width="196" height="34" rx="8" fill="#fff" stroke="#d5deea"/>
<text class="m" x="572" y="96">1</text><text class="s" x="572" y="109">intercept</text>
<rect x="560" y="123" width="196" height="34" rx="8" fill="#fff" stroke="#d5deea"/>
<text class="m" x="572" y="138">K<tspan baseline-shift="sub" font-size="9">ij</tspan></text><text class="s" x="572" y="151">genetic relatedness (GRM)</text>
<rect x="560" y="165" width="196" height="34" rx="8" fill="#fff" stroke="#d5deea"/>
<text class="m" x="572" y="180">1(i = j)</text><text class="s" x="572" y="193">same donor</text>
<rect x="560" y="207" width="196" height="34" rx="8" fill="#fff" stroke="#d5deea"/>
<text class="m" x="572" y="222">c<tspan baseline-shift="sub" font-size="9">im</tspan> c<tspan baseline-shift="sub" font-size="9">jn</tspan></text><text class="s" x="572" y="235">cell state, covariates</text>
<rect x="560" y="249" width="196" height="34" rx="8" fill="#e6f4f4" stroke="#009999"/>
<text class="m" x="572" y="264">K<tspan baseline-shift="sub" font-size="9">ij</tspan> c<tspan baseline-shift="sub" font-size="9">im</tspan> c<tspan baseline-shift="sub" font-size="9">jn</tspan></text><text class="s" x="572" y="277">G×cell state</text>
<text class="s" x="560" y="314">least squares on donor-level sums</text>
<path d="M766 170 L786 170" stroke="#0074d9" stroke-width="2" marker-end="url(#ar)"/>
<rect x="800" y="221.6" width="44" height="68.4" fill="#dbe3ee"/>
<text class="l" x="854" y="259.8">noise</text>
<rect x="800" y="203.4" width="44" height="18.2" fill="#9aa7b8"/>
<text class="l" x="854" y="216.5">covariate</text>
<rect x="800" y="176.0" width="44" height="27.4" fill="#0b2a5c"/>
<text class="l" x="854" y="193.7">donor I</text>
<rect x="800" y="125.8" width="44" height="50.2" fill="#009999"/>
<text class="l" x="854" y="154.9">G×cell state</text>
<rect x="800" y="62.0" width="44" height="63.8" fill="#0074d9"/>
<text class="l" x="854" y="97.9">G</text>
<rect x="800" y="62" width="44" height="228" rx="6" fill="none" stroke="#fff" stroke-width="3"/>
<text class="s" x="800" y="312">expression variance</text>
</svg>
</div>


<div class="sc-cards">
<a class="sc-card" href="articles/tutorial-cell-state.html">
<span class="sc-tag sc-tag-blue">schint_cell()</span>
<h3>G×Cell-State</h3>
<p>Main genetic and genotype × cell-state heritability from single-cell expression of one cell type.</p>
</a>
<a class="sc-card" href="articles/tutorial-cell-type.html">
<span class="sc-tag sc-tag-teal">schint_sample()</span>
<h3>G×Cell-Type</h3>
<p>Cell-type-shared and cell-type-specific genetic variance from donor × cell-type pseudobulk.</p>
</a>
<a class="sc-card" href="articles/tutorial-perturbation.html">
<span class="sc-tag sc-tag-navy">schint_pert()</span>
<h3>Perturbation × context</h3>
<p>Perturbation and perturbation × cell-state / cell-type variance in perturb-seq screens.</p>
</a>
</div>





## 1. Introduction

**scHINT** (single-cell gene expression **H**eritability with gene-by-cell-context
**INT**eractions) quantifies how much of the variation in a gene's expression is
explained by genetics, and how much of that genetic effect **depends on the cellular
context**. Most expression-QTL and heritability analyses use pseudobulk profiles, which
average over cells and so cannot see genetic effects that change along a continuous cell
state. scHINT works directly on single-cell data.

**The main model is the cell-level G×Cell-State heritability model**
(`schint_cell()`): for each gene, the variance of single-cell expression is partitioned into

* a **main genetic** component (G), shared by all cells,
* a **G×Cell-State** component: genetic effects that vary continuously along cell-state axes
  such as principal components or pseudotime,
* a donor (individual) component, cell-state and covariate components, and cell-level noise.

Estimation is a Haseman-Elston (HE) regression on pairs of cells that is assembled from
donor-level sufficient statistics, so it scales to millions of cells, and standard errors
come from a donor-block jackknife at almost no extra cost.

The same framework extends to other contexts:

| model | function | what it estimates |
|---|---|---|
| **G×Cell-State** (cell level) | `schint_cell()` | main genetic and genotype × cell-state heritability from single-cell data |
| **G×Cell-Type** (pseudobulk) | `schint_sample()` | cell-type-shared and cell-type-specific genetic variance from donor × cell-type expression |
| **Perturbation × context** | `schint_pert()` | perturbation and perturbation × cell-state / cell-type variance in perturb-seq screens |

## 2. Models

Full derivations, including the
sufficient-statistics algorithm and the jackknife, are on the
[Model page](https://LidaWangPSU.github.io/scHINT/articles/model.html).

### 2.1 G×Cell-State heritability (cell level)

For a gene, let $y_{im}$ be the normalized expression in cell $m$ ($m=1,\dots,M_i$) of
individual $i$ ($i=1,\dots,n$), $g_i\in\mathbb{R}^p$ the standardized cis-genotypes (variants
within ±500 kb of the gene), $c_{im}\in\mathbb{R}^K$ the cell-state features (e.g. the leading
within-cell-type principal components) and $o_{im}$ further covariates:

$$
y_{im}=g_i^\top\beta+(g_i\otimes c_{im})^\top\gamma+c_{im}^\top\alpha+o_{im}^\top\eta+e_i+\varepsilon_{im},
$$

$$
\beta\sim N\!\left(0,\sigma_G^2p^{-1}I_p\right),\quad
\gamma_k\sim N\!\left(0,\sigma_{G\times C,k}^2p^{-1}I_p\right),\quad
e_i\sim N(0,\sigma_I^2),\quad
\varepsilon_{im}\sim N(0,\sigma_\varepsilon^2).
$$

$\sigma_G^2$ is the contribution of main cis-genetic effects and $\sigma_{G\times C,k}^2$ that of
genetic effects varying along cell-state axis $k$. With $K^G_{ij}=g_i^\top g_j/p$ the cis-GRM, the expected
product of two different cells is

$$
\mathbb{E}(y_{im}y_{jn})=\sigma_G^2K^G_{ij}+\sum_{k=1}^{K_{\rm int}}\sigma_{G\times C,k}^2K^G_{ij}c_{imk}c_{jnk}
+\sigma_I^2\mathbf 1(i=j)+\sum_{k=1}^{K}\sigma_{C,k}^2c_{imk}c_{jnk}+\sum_{l=1}^{L}\sigma_{O,l}^2o_{iml}o_{jnl}.
$$

scHINT regresses $y_{im}y_{jn}$ on these kernels. The baseline model

$$
y_{im}y_{jn}\sim 1+K^G_{ij}+\mathbf 1(i=j)+\sum_{k=1}^{K}c_{imk}c_{jnk}+\sum_{l=1}^{L}o_{iml}o_{jnl}
$$

estimates the main genetic component, and the interaction model adds
$\sum_{k=1}^{K_{\rm int}}K^G_{ij}c_{imk}c_{jnk}$; the coefficients of $K^G_{ij}$ and
$K^G_{ij}c_{imk}c_{jnk}$ are $\hat\sigma_G^2$ and $\hat\sigma_{G\times C,k}^2$, and
$\hat\sigma_{G\times C}^2=\sum_k\hat\sigma_{G\times C,k}^2$.
To make the estimates comparable to population-level heritability, which is not diluted by
cell-level noise, heritability can be defined with the cell-level residual removed from the
denominator:

$$
h_G^2=\frac{\sigma_G^2}{\sigma_G^2+\sigma_{G\times C}^2+\sigma_C^2+\sigma_o^2+\sigma_I^2},\qquad
h_{G\times C}^2=\frac{\sigma_{G\times C}^2}{\sigma_G^2+\sigma_{G\times C}^2+\sigma_C^2+\sigma_o^2+\sigma_I^2}.
$$

A list of GRMs (e.g. MAF bins) gives one $\sigma_G^2$ and one $\sigma_{G\times C}^2$ per GRM.

### 2.2 G×Cell-Type heritability (pseudobulk)

For the pseudobulk expression $y_{it}$ of individual $i$ in cell type $t=1,\dots,T$:

$$
y_{it}=g_i^\top\beta+g_i^\top\gamma_t+\alpha_t+o_{it}^\top\eta+e_i+\varepsilon_{it},
$$

$$
\beta\sim N\!\left(0,\sigma_G^2p^{-1}I_p\right),\quad
\gamma_t\sim N\!\left(0,\sigma_{G\times CT}^2p^{-1}I_p\right),\quad
\alpha_t\sim N(0,\sigma_T^2),\quad e_i\sim N(0,\sigma_I^2),
$$

where $\beta$ is the genetic effect shared across cell types, $\gamma_t$ the cell-type-specific
deviation, $\alpha_t$ the cell-type mean and $e_i$ a donor effect shared across cell types.
For two observations $(i,s)\neq(j,t)$

$$
\mathbb{E}(y_{is}y_{jt})=\sigma_G^2K^G_{ij}+\sigma_{G\times CT}^2K^G_{ij}\mathbf 1(s=t)
+\sigma_{CT}^2\mathbf 1(s=t)+\sigma_I^2\mathbf 1(i=j),
$$

and the HE regression uses the kernels $K^G_{ij}$, $K^G_{ij}\mathbf 1(s=t)$, $\mathbf 1(s=t)$ and
$\mathbf 1(i=j)$ (plus covariate kernels). $\sigma_G^2$ is genetic variance shared across cell types and
$\sigma_{G\times CT}^2$ is cell-type-specific genetic variance.

### 2.3 Perturbation × context heritability

For a gene, let $p_m$ be the perturbation of cell $m$ and $c_m\in\mathbb{R}^K$ its standardized cell-state
features:

$$
y_m=b_{p_m}+\sum_{k=1}^{K}c_{mk}\alpha_k+\sum_{k=1}^{K}c_{mk}u_{p_mk}+\varepsilon_m,
$$

$$
b_{p}\sim N(0,\sigma_P^2),\quad \alpha_k\sim N(0,\sigma_{{\rm PC},k}^2),\quad
u_{pk}\sim N(0,\sigma_{P\times{\rm PC},k}^2),\quad \varepsilon_m\sim N(0,\sigma_\varepsilon^2).
$$

$\sigma_P^2$ is perturbation-associated variation shared across cell states and
$\sigma_{P\times{\rm PC},k}^2$ the variance of the perturbation response along axis $k$. For two
cells $m\neq n$ the HE regression is

$$
y_my_n\sim 1+\mathbf 1(p_m=p_n)+\sum_{k=1}^{K}c_{mk}c_{nk}+\mathbf 1(p_m=p_n)\sum_{k=1}^{K}c_{mk}c_{nk}.
$$

With a categorical context (cell type) $c_{mk}c_{nk}$ is replaced by $\mathbf 1(t_m=t_n)$. Because genes measured
in the same cells share their pair predictors, the design cross-products are computed once for all genes.

## 3. Install

```r
install.packages("remotes")
remotes::install_github("LidaWangPSU/scHINT")
```

scHINT needs only base R. (`BEDMatrix` is optional, to read PLINK files with `read_plink()`.)
The examples below use the simulated data shipped with the package:


``` r
data(schint_example)        # genotypes + cell-level and sample-level expression
ex <- schint_example
```

## 4. Usage

### 4.1 G×Cell-State heritability: `schint_cell()`

**Input**

* `data` -- one row per cell (of one cell type) with expression, donor ID, cell-state
  variables and covariates.
* `y`, `id` -- names of the expression and donor-ID columns.
* `geno` -- donors × cis-SNPs dosage matrix (row names = donor IDs); one GRM is built from
  all SNPs. **Or** `grm` -- a pre-computed GRM, or a list of GRMs for a multi-GRM model.
* `context` -- cell-state column(s) that interact with genetics (numeric = continuous,
  factor = categorical).
* `covariates` -- other covariates (main effects only). `jackknife = TRUE` for standard errors.


``` r
head(ex$cell[, c("donor", "state", "pc1", "pc2", "age", "sex", "batch", "GENE_A")], 3)
#>   donor       state        pc1        pc2 age sex batch     GENE_A
#> 1  D001 -0.34413687 -0.6334822 0.02127604  63   F    B1 -1.1423329
#> 2  D001  0.04947898 -0.7615181 0.84821092  63   F    B1 -0.3653666
#> 3  D001 -0.84762422  1.2939331 0.24429519  63   F    B1  0.7018750
dim(ex$geno)
#> [1] 400 250
```

**Example**


``` r
fit <- schint_cell(ex$cell, y = "GENE_A", id = "donor", geno = ex$geno,
                   context = "state",
                   covariates = c("pc1", "pc2", "age", "sex", "batch"),
                   jackknife = TRUE, n_blocks = 50)
fit
#> scHINT model: 12000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2      se   se_h2
#>                G  0.20500 0.4140 0.04700 0.08190
#>              GxC  0.16300 0.3290 0.03890 0.04980
#>          G_total  0.36900 0.7430 0.06680 0.07080
#>                I  0.07670 0.1550 0.03580 0.07430
#>          context  0.00704 0.0142 0.00331 0.00617
#>        covariate  0.04360 0.0879 0.01130 0.02600
#>  total_explained  0.49600     NA 0.06100      NA
#> 
#> Standard errors: jackknife over 50 donor blocks
```

**Output** -- `fit$summary` (grouped components) and `fit$coefficients` (one row per
regression term): `G` main genetic variance, `GxC` G×Cell-State variance, `G_total = G + GxC`,
`I` donor, `context` and `covariate` components; `h2` population-level heritability (residual removed from the denominator);
jackknife `se*` columns. Several GRMs and several context variables add per-GRM and
per-variable rows.


``` r
# several GRMs: e.g. cis SNPs split into three MAF bins, and two cell-state axes
grms <- lapply(split_snps(ex$geno, 3, by = "maf"), make_grm)
schint_cell(ex$cell, "GENE_A", "donor", grm = grms,
            context = c("state", "subtype"), covariates = c("pc1", "pc2"))$summary[, 1:3]
#>               term    variance    fraction
#> 1                G  0.21304077  0.21304077
#> 2                1  0.06190295  0.06190295
#> 3                2  0.05429398  0.05429398
#> 4                3  0.09684385  0.09684385
#> 5              GxC  0.15004511  0.15004511
#> 6              1:C  0.05324709  0.05324709
#> 7              2:C  0.04827433  0.04827433
#> 8              3:C  0.04852369  0.04852369
#> 9          G_total  0.36308588  0.36308588
#> 10       GxC:state  0.16876151  0.16876151
#> 11     GxC:subtype -0.01871639 -0.01871639
#> 12               I  0.09303695  0.09303695
#> 13         context  0.02630493  0.02630493
#> 14       covariate  0.02027540  0.02027540
#> 15 total_explained  0.50270316  0.50270316
```

### 4.2 G×Cell-Type heritability (pseudobulk): `schint_sample()`

**Input** -- `data` with one row per donor × cell type (expression, donor ID, cell-type
column, optional covariates); the same genetic input as above; `celltype` names the
cell-type column.


``` r
head(ex$sample[, c("donor", "celltype", "age", "sex", "GENE_A")], 3)
#>   donor celltype age sex     GENE_A
#> 1  D001      CT1  63   F -1.2743896
#> 2  D001      CT2  63   F -0.4401739
#> 3  D001      CT3  63   F -2.2811747
```

**Example**


``` r
fs <- schint_sample(ex$sample, y = "GENE_A", id = "donor", celltype = "celltype",
                    geno = ex$geno, covariates = c("age", "sex"), jackknife = TRUE)
fs
#> scHINT model: 2000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2     se  se_h2
#>                G   0.1020 0.1080 0.0237 0.0248
#>              GxC   0.0566 0.0601 0.0206 0.0204
#>          G_total   0.1580 0.1680 0.0295 0.0279
#>                I   0.0228 0.0243 0.0199 0.0211
#>          context   0.7320 0.7780 0.0216 0.0215
#>        covariate   0.0278 0.0295 0.0103 0.0106
#>  total_explained   0.9410     NA 0.0362     NA
#> 
#> Standard errors: jackknife over 400 donor blocks
```

**Output** -- same structure as above: `G` is genetic variance shared across cell types,
`GxC` the cell-type-specific genetic variance, `I` the donor component shared across cell
types, `context` the cell-type mean variance. `gxc = FALSE` fits shared genetics only;
`cat_mode = "per_level"` estimates a genetic variance for each cell type.

### 4.3 Perturbation × context heritability: `schint_pert()`

**Input** -- `data` with one row per cell, one or several expression columns (`y`), the
`perturb` column, `context` (cell state and/or cell type), `covariates`, `control` labels to drop
(e.g. `"NT"`) and optionally `perturb_group`.


``` r
data(schint_pert_example)
head(schint_pert_example[, c("perturb", "perturb_group", "celltype", "state", "GENE_A")], 3)
#>   perturb perturb_group celltype     state     GENE_A
#> 1  PERT01        group1      CT2 1.1188078 0.07107567
#> 2  PERT01        group1      CT3 1.1119961 0.14960499
#> 3  PERT01        group1      CT1 0.4430754 0.38410583
```

**Example**


``` r
fp <- schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B"), perturb = "perturb",
                  context = c("state", "celltype"), covariates = c("pc1", "pc2"),
                  control = "NT", jackknife = TRUE)
fp
#> scHINT perturbation model: 3470 cells, 60 perturbations, 2 genes
#> 
#>    gene            term variance      se
#>  GENE_A               P  0.02720 0.00866
#>  GENE_A             PxC  0.06180 0.01930
#>  GENE_A         P_total  0.08900 0.02080
#>  GENE_A       PxC:state  0.05070 0.01460
#>  GENE_A    PxC:celltype  0.01110 0.01130
#>  GENE_A         context  0.09800 0.01260
#>  GENE_A       covariate  0.04050 0.00617
#>  GENE_A total_explained  0.22700 0.02510
#>  GENE_B               P  0.04350 0.01510
#>  GENE_B             PxC  0.02450 0.01220
#>  GENE_B         P_total  0.06800 0.01880
#>  GENE_B       PxC:state -0.00277 0.00293
#>  GENE_B    PxC:celltype  0.02730 0.01260
#>  GENE_B         context  0.06640 0.01250
#>  GENE_B       covariate  0.05610 0.00832
#>  GENE_B total_explained  0.19100 0.02470
#> 
#> Standard errors: jackknife over 60 perturbation blocks
```

**Output** -- long-format `fp$summary` / `fp$coefficients` with one block per gene: `P`
perturbation variance shared across states, `PxC` perturbation × context variance (split
into `PxC:state` and `PxC:celltype`), `P_total`, and per-group rows when `perturb_group` is used.

### More

[Users manual](https://LidaWangPSU.github.io/scHINT/articles/users-manual.html) ·
[Input formats and outputs](https://LidaWangPSU.github.io/scHINT/articles/inputs-outputs.html) ·
[Tutorials](https://LidaWangPSU.github.io/scHINT/articles/tutorial-cell-state.html) ·
[Troubleshooting](https://LidaWangPSU.github.io/scHINT/articles/troubleshooting.html)

**Citation.** Wang L, Fadil C, Wang L, Van Dyke K, Gusev A. *Single-cell gene expression heritability
informed by gene × cell-context interactions.* (manuscript)
